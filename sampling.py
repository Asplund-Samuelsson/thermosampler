#!/usr/bin/env python3

# Import libraries
import numpy as np
import pandas as pd
import collections
import argparse
import re
import sys, os

# Import functions from MDF script
from mdf import parse_equation, read_reactions, read_reaction_drGs, remove_comments
from mdf import read_constraints, sWrite, sError
from mdf import read_ratio_constraints

# Define function for random sampling of concentrations
def random_c(c_lim):
    sample = np.array([np.random.random() for n in range(0, c_lim.shape[0])])
    return sample * (c_lim[:,1] - c_lim[:,0]) + c_lim[:,0]

# Define function for checking if set is thermodynamically feasible
def df_ok(c, g, RT, S, mdf=0):
    # Calculate delta G prime
    df = -(g + RT * np.sum(np.transpose(S) * c, 1))
    # Check if all driving forces higher than the minimum driving force
    return sum(df > mdf) == df.shape[0]

# Define function for checking if set has acceptable ratios
def ratios_ok(c, ratio_lim, ratio_mat):
    ratios = np.sum(ratio_mat.T * c, 1).reshape([ratio_lim.shape[0], 1])
    min = np.sum(np.subtract(ratios, ratio_lim) > 0, 0)[0] == ratios.shape[0]
    max = np.sum(np.subtract(ratios, ratio_lim) < 0, 0)[1] == ratios.shape[0]
    return min and max

# Define function for checking that sum of concentrations is not too high
def sum_ok(c, max_tot_c):
    return np.sum(np.exp(c)) <= max_tot_c

# Define function for checking custom metabolite group concentration sums are ok
def sums_ok(c, c_sums):
    # If there are no concentration sum constraints, return True
    if c_sums is None:
        return True
    # Check sum of each group
    for i in c_sums[0].keys():
        if np.sum(np.exp(c[list(c_sums[0][i])])) > c_sums[1][i]:
            return False
    # If all groups passed, return True
    return True

# Define function that checks concentrations are within limits
def limits_ok(c, c_lim):
    c_l = c.reshape([c.shape[0],1])
    min = np.sum(np.subtract(c_l, c_lim) >= 0, 0)[0] == c.shape[0]
    max = np.sum(np.subtract(c_l, c_lim) <= 0, 0)[1] == c.shape[0]
    return min and max

# Define function for checking feasibility, ratios, sum, and limits in one go
def is_feasible(
    c, g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim, c_sums, mdf=0
    ):
    # Run all feasibility tests
    ok = [
        # Driving forces above minimum
        df_ok(c, g, RT, S, mdf),
        # Concentration ratios within bounds
        ratios_ok(c, ratio_lim, ratio_mat),
        # Concentration sum below max
        sum_ok(c, max_tot_c),
        # Concentrations within bounds
        limits_ok(c, c_lim),
        # Concentration sum groups below individual max values
        sums_ok(c, c_sums)
    ]
    # Check that all feasibility tests are True
    return sum(ok) == len(ok)

# Define function for selecting a random direction
def random_direction(c_lim):
    # Create a random vector of the same length as c
    direction = np.array([np.random.random() for n in range(0, c_lim.shape[0])])
    # Subtract 0.5 to introduce negative directions
    direction = direction - 0.5
    # Set fixed concentration direction to zero
    direction[c_lim[:,1] - c_lim[:,0] == 0] = 0
    # Normalize length of direction vector
    normalized_direction = direction / np.linalg.norm(direction)
    return normalized_direction

# Define function to generate one feasible metabolite concentration set
def generate_feasible_c(
        g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim, mdf=0, c_sums=None
    ):
    c = random_c(c_lim) # Initialize c
    while not is_feasible(
        c, g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim, c_sums, mdf
    ):
        c = random_c(c_lim) # Generate new c until feasible
    return c

# Modify direction in order to get unstuck from concentration limits
def unstick_direction(c, direction, c_lim):
    # Determine what metabolites are stuck at limits
    stuck = c.reshape((c.size,1)) == c_lim
    # Determine current signs of direction vector
    dirsign = np.sign(direction)
    # Pick a random sign for metabolites stuck at max
    max_sign = np.random.choice([-1,1], 1)
    # All directions for metabolites stuck at max must be the same sign
    dirsign[stuck[:,1] * dirsign != 0] = max_sign
    # All directions for metabolites stuck at min must be the opposite sign
    dirsign[stuck[:,0] * dirsign != 0] = -max_sign
    # Determine the directions that must change sign
    change_sign = dirsign != np.sign(direction)
    # Change the sign of directions that must change sign
    direction[change_sign] = direction[change_sign] * -1
    # Return the compatibility-modified "unstuck" direction vector
    return direction

# Determine minimum and maximum possible theta given concentration limits
def calculate_theta_hard_limit(c, direction, c_lim):
    # Find smallest fraction of direction that hits limit if added
    theta_max = np.vstack([
        (c_lim[:,1] - c)[direction != 0] / direction[direction != 0],
        (c_lim[:,0] - c)[direction != 0] / direction[direction != 0]
    ])
    theta_max = np.max(theta_max, 0)
    theta_max = min(theta_max[theta_max >= 0])
    # Find smallest fraction of direction that hits limit if subtracted
    theta_min = np.vstack([
        (c - c_lim[:,1])[direction != 0] / direction[direction != 0],
        (c - c_lim[:,0])[direction != 0] / direction[direction != 0]
    ])
    theta_min = np.max(theta_min, 0)
    theta_min = -min(theta_min[theta_min >= 0])
    return (theta_min, theta_max)

# Define function for determining minimum and maximum step length (theta)
def theta_range(
    c, g, RT, S, ratio_lim, ratio_mat, max_tot_c,
    c_lim, direction, precision=1e-3, mdf=0, c_sums=None
    ):
    # Define function for honing in on a theta limit
    def hone_theta(theta_outer, theta_inner=0):
        if is_feasible(
            c + theta_outer * direction, g, RT, S, ratio_lim, ratio_mat,
            max_tot_c, c_lim, c_sums, mdf
            ):
            # If the outer theta is feasible, accept that solution
            theta_inner = theta_outer
        else:
            while abs(theta_outer - theta_inner) > precision:
                # Calculate a theta between outer and inner limits
                theta_cur = (theta_outer + theta_inner) / 2
                if is_feasible(
                    c + theta_cur * direction, g, RT, S, ratio_lim, ratio_mat,
                    max_tot_c, c_lim, c_sums, mdf
                    ):
                    # Move outwards, set inner limit to current theta
                    theta_inner = theta_cur
                else:
                    # Move inwards, set outer limit to current theta
                    theta_outer = theta_cur
        # Return inner theta
        return theta_inner
    # Get hard limits on theta from concentrations
    theta_lim = calculate_theta_hard_limit(c, direction, c_lim)
    # Hone in on upper theta
    theta_upper = hone_theta(theta_lim[1])
    # Hone in on lower theta
    theta_lower = hone_theta(theta_lim[0])
    # Return results
    return [theta_lower, theta_upper]

# Define function for performing hit-and-run sampling within the solution space
def hit_and_run(
    c, g, RT, S, ratio_lim, ratio_mat, max_tot_c,
    c_lim, n_samples, precision=1e-3, mdf=0, n=1, c_sums=None
    ):
    # Set up concentration storage list
    fMCSs = [c]
    # Perform n steps
    for i in range(0, n_samples):
        # Generate random direction
        direction = random_direction(c_lim)
        # If it is the first step, unstick the direction
        if i == 0:
            direction = unstick_direction(c, direction, c_lim)
        # Determine minimum and maximum step length
        theta = theta_range(
            c, g, RT, S, ratio_lim, ratio_mat, max_tot_c,
            c_lim, direction, precision, mdf, c_sums
        )
        # Perform a random sampling of the step length
        theta = theta[0] + np.random.random() * (theta[1] - theta[0])
        # Perform step
        c = c + theta * direction
        # Ensure feasibility
        if not is_feasible(
                c, g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim, c_sums, mdf
            ):
            print("Warning: Infeasible point reached.")
            break
        # Store concentration if the step is correct
        if (i + 1) % n == 0:
            fMCSs.append(c)
    # Return list of concentrations
    return fMCSs

def make_ratio_mat(ratio_constraints, S_pd):
    # Prepare ratio matrix of correct shape with only zeros
    rat_mat = np.zeros((S_pd.shape[0], ratio_constraints.shape[0]))
    # Iterate over ratio indices
    for rat_i in ratio_constraints.index:
        # Add the numerator
        rat_mat[
            np.where(S_pd.index == ratio_constraints['cpd_id_num'][rat_i]),
            rat_i
        ] = 1
        # Add the denominator
        rat_mat[
            np.where(S_pd.index == ratio_constraints['cpd_id_den'][rat_i]),
            rat_i
        ] = -1
    # Return finished ratio matrix
    return rat_mat

def read_concentrations(c_text):
    if "," in c_text:
        # Parse MDF table (metabolite columns)
        c_text = [x.split(',') for x in c_text.strip().split("\n")]
        i = np.where([x.startswith('c_') for x in c_text[0]])[0]
        concentrations = pd.DataFrame(dict(zip(
            ['Conc' + str(x) for x in range(1,len(c_text))],
            [[float(c_text[r][c]) for c in i] for r in range(1,len(c_text))]
        )))
        concentrations.index = [c_text[0][x].strip('c_') for x in i]
    else:
        # Parse tab-delimited table (metabolite rows)
        c_text = [x.split('\t') for x in c_text.strip().split("\n")]
        concentrations = pd.DataFrame(dict(zip(
            ['Conc' + str(x) for x in range(1,len(c_text[0]))],
            map(list, zip(*[
                [float(y) for y in c_text[x][1:]] for x in range(len(c_text))
            ]))
        )))
        concentrations.index = [c_text[x][0] for x in range(len(c_text))]
    return concentrations

def read_sums(sums_text, S_pd):
    # Parse text
    sums_list = [remove_comments(x).split("\t") for x in sums_text.strip().split("\n")]
    # Remove empty entries (from blank lines)
    sums_list = [x for x in sums_list if len(x) > 1]
    # Make into data frame
    sums_df = pd.DataFrame(sums_list)
    # Find row index of each compound in stoichiometric matrix
    c_sums = sums_df.copy()
    c_sums.loc[:,0] = [S_pd.index.tolist().index(x) for x in sums_df.loc[:,0]]
    # Rename columns
    c_sums = c_sums.rename(columns={0:'cpd_index',1:'cpd_group',2:'group_sum'})
    # Convert group and sum to int and float
    c_sums.loc[:,'cpd_group'] = pd.to_numeric(c_sums['cpd_group'])
    c_sums.loc[:,'group_sum'] = pd.to_numeric(c_sums['group_sum'])
    # Error and terminate if more than one sum per group
    c_sums_grouped = c_sums.drop('cpd_index', axis=1).groupby('cpd_group')
    c_sums_check = c_sums_grouped.agg(lambda x: len(set(x)))['group_sum']
    for x in c_sums_check.tolist():
        if x > 1:
            sError("\nError: More than one concentration sum per group.\n")
            exit()
    # Modify concentration sums table to be more lean
    c_groups = {}
    for i in range(c_sums.shape[0]):
        cpd_index = c_sums['cpd_index'].tolist()[i]
        try:
            c_groups[c_sums['cpd_group'].tolist()[i]].add(cpd_index)
        except KeyError:
            c_groups[c_sums['cpd_group'].tolist()[i]] = {cpd_index}
    c_sums = (
        c_groups,
        dict(zip(c_sums['cpd_group'].tolist(), c_sums['group_sum'].tolist()))
    )
    return c_sums

# Main code block
def main(
        reaction_file, std_drG_file, outfile_name, cons_file, ratio_cons_file,
        max_tot_c, n_samples, n_starts=1, c_file = None,
        T=298.15, R=8.31e-3, proton_name='C00080', precision=1e-3, mdf=0, n=1,
        conc_sums=None
    ):

    sWrite("\nLoading data...")

    RT = R*T

    # Load stoichiometric matrix
    S_pd = read_reactions(open(reaction_file, 'r').read(), proton_name)
    S = S_pd.to_numpy()

    # Load standard reaction Gibbs energies
    std_drGs = read_reaction_drGs(open(std_drG_file, 'r').read())
    std_drGs.index = std_drGs['rxn_id']
    std_drGs = std_drGs.reindex(S_pd.columns.values)
    g = std_drGs['drG'].to_numpy()

    # Load concentration constraints
    constraints = read_constraints(open(cons_file, 'r').read())
    constraints.index = constraints['cpd_id']
    constraints = constraints.reindex(S_pd.index.values)
    c_lim = np.log(constraints[['x_min', 'x_max']].to_numpy())

    # Load concentration ratio constraints
    if ratio_cons_file:
        ratio_cons_text = open(ratio_cons_file, 'r').read()
    else:
        ratio_cons_text = ''
    ratio_constraints = read_ratio_constraints(ratio_cons_text)
    ratio_mat = make_ratio_mat(ratio_constraints, S_pd)
    ratio_lim = np.log(ratio_constraints[['ratio','ratio_upper']].to_numpy())

    # Load concentrations
    if not c_file:
        c_loaded = False
    else:
        c_loaded = True
        c_pd = read_concentrations(open(c_file, 'r').read())
        c_pd = c_pd.reindex(S_pd.index.values)

    # Load concentration sums
    if not conc_sums:
        c_sums = None
    else:
        c_sums = read_sums(open(conc_sums, 'r').read(), S_pd)

    sWrite(" Done.\n")

    # Make sure we can write to outfile before the heavy lifting
    try:
        outfile = open(outfile_name, 'w')
    except IOError as e:
        print('Unable to open outfile for writing:', str(e))
        sys.exit()

    sWrite("Performing hit-and-run sampling...")
    if c_loaded:
        fMCSs = [
            hit_and_run(
                np.log(c_pd.iloc[:,i].to_numpy()),
                g, RT, S, ratio_lim, ratio_mat, max_tot_c,
                c_lim, n_samples, precision, mdf, n, c_sums
            ) \
            for i in range(c_pd.shape[1])
        ]
    else:
        fMCSs = [
            hit_and_run(
                generate_feasible_c(g,RT,S,ratio_lim,ratio_mat,max_tot_c,c_lim),
                g, RT, S, ratio_lim, ratio_mat, max_tot_c,
                c_lim, n_samples, precision, mdf, n, c_sums
            )\
            for i in range(n_starts)
        ]
    sWrite(" Done.\n")

    # Save data to outfile
    header = 'Run\tfMCS\t' + "\t".join(S_pd.index.values) + "\n"
    junk = outfile.write(header)
    r = 0
    for run in fMCSs:
        f = 0
        for fMCS in run:
            fMCS_string = "\t".join([str(x) for x in fMCS])
            output = "\t".join([str(r), str(f)]) + "\t" + fMCS_string + "\n"
            junk = outfile.write(output)
            f += n
        r +=1

    outfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--reactions', type=str, required=True,
        help='Load reactions.'
        )
    parser.add_argument(
        '--std_drG', type=str, required=True,
        help='Load standard reaction Gibbs energies.'
        )
    parser.add_argument(
        '--outfile', type=str, required=True,
        help='Write fMCSs table in tab-delimited format.'
        )
    parser.add_argument(
        '--constraints', type=str, required=True,
        help='Load concentration bound constraints.'
        )
    parser.add_argument(
        '--ratios', type=str, default=None,
        help='Load concentration ratio constraints.'
        )
    parser.add_argument(
        '--concs', type=str, default=None,
        help='Load start concentrations in tab or MDF csv format (optional).'
        )
    parser.add_argument(
        '--steps', type=int, default=1000,
        help="Number of hit-and-run steps."
        )
    parser.add_argument(
        '--starts', type=int, default=1,
        help="Number of starts if no starting concentrations (default 1)."
        )
    parser.add_argument(
        '-T', type=float, default=298.15,
        help='Temperature (K).'
        )
    parser.add_argument(
        '-R', type=float, default=8.31e-3,
        help='Universal gas constant (kJ/(mol*K)).'
        )
    parser.add_argument(
        '--proton_name', default='C00080',
        help='Name used to identify protons.'
        )
    parser.add_argument(
        '--max_conc', type=float, default=1.05,
        help='Maximum total concentration (M).'
        )
    parser.add_argument(
        '--precision', type=float, default=1e-3,
        help='Precision of step size calculation (default 1e-3).'
        )
    parser.add_argument(
        '--mdf', type=float, default=0,
        help='Minimum driving force (default 0).'
        )
    parser.add_argument(
        '-n', type=int, default=1,
        help="Save every nth step of the random walk (default 1)."
    )
    parser.add_argument(
        '--conc_sums', default=None,
        help='Tab-delimited compound name, group, and group sum (default None).'
    )
    args = parser.parse_args()
    main(
        args.reactions, args.std_drG, args.outfile, args.constraints,
        args.ratios, args.max_conc, args.steps, args.starts,
        args.concs, args.T, args.R, args.proton_name,
        args.precision, args.mdf, args.n, args.conc_sums
    )
