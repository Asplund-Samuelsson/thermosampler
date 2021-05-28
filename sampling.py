# Import libraries
import numpy as np

# Stoichiometric matrix (TCA cycle); columns are reactions
# R00351	C00024 + C00001 + C00036 <=> C00158 + C00010
# R01325	C00158 <=> C00417 + C00001
# R01900    C00001 + C00417 <=> C00311
# R00709	C00311 + C00003 <=> C00026 + C00011 + C00004 + C00080
# R08549	C00026 + C00010 + C00003 <=> C00091 + C00011 + C00004 + C00080
# R00405	C00008 + C00009 + C00091 <=> C00002 + C00042 + C00010
# M00148	C00399 + C00042 <=> C00390 + C00122
# R01082	C00122 + C00001 <=> C00149
# R00342	C00149 + C00003 <=> C00036 + C00004 + C00080
S = np.array([
    # R00351 R01325 R01900 R00709 R08549 R00405 M00148 R01082 R00342
    [     -1,     1,    -1,     0,     0,     0,     0,    -1,     0],  # C00001
    [      0,     0,     0,     0,     0,     1,     0,     0,     0],  # C00002
    [      0,     0,     0,    -1,    -1,     0,     0,     0,    -1],  # C00003
    [      0,     0,     0,     1,     1,     0,     0,     0,     1],  # C00004
    [      0,     0,     0,     0,     0,    -1,     0,     0,     0],  # C00008
    [      0,     0,     0,     0,     0,    -1,     0,     0,     0],  # C00009
    [      1,     0,     0,     0,    -1,     1,     0,     0,     0],  # C00010
    [      0,     0,     0,     1,     1,     0,     0,     0,     0],  # C00011
    [     -1,     0,     0,     0,     0,     0,     0,     0,     0],  # C00024
    [      0,     0,     0,     1,    -1,     0,     0,     0,     0],  # C00026
    [     -1,     0,     0,     0,     0,     0,     0,     0,     1],  # C00036
    [      0,     0,     0,     0,     0,     1,    -1,     0,     0],  # C00042
#   [      0,     0,     0,     0,     0,     0,     0,     0,     0],  # C00080
    [      0,     0,     0,     0,     1,    -1,     0,     0,     0],  # C00091
    [      0,     0,     0,     0,     0,     0,     1,    -1,     0],  # C00122
    [      0,     0,     0,     0,     0,     0,     0,     1,    -1],  # C00149
    [      1,    -1,     0,     0,     0,     0,     0,     0,     0],  # C00158
    [      0,     1,    -1,     0,     0,     0,     0,     0,     0],  # C00417
    [      0,     0,     1,    -1,     0,     0,     0,     0,     0],  # C00311
    [      0,     0,     0,     0,     0,     0,     1,     0,     0],  # C00390
    [      0,     0,     0,     0,     0,     0,    -1,     0,     0]   # C00399
])

# Delta G values @ pH 7.6 (TCA cycle)
# R00351	-38.1
# R01325	8.3
# R01900    -0.7
# R00709	5.4
# R08549	-27.2
# R00405	1.1
# M00148	-21.6
# R01082	-3.4
# R00342	27.0
g = np.array([-38.1, 8.3, -0.7, 5.4, -27.2, 1.1, -21.6, -3.4, 27.0])

# RT is a constant
T=298.15
R=8.31e-3
RT = R*T

# Constrain concentrations
# Use natural log
c_lim = np.log(np.array([
    [           1,           1 ], # C00001
    [ 0.000210837,    0.040093 ], # C00002
    [ 0.000386754,  0.00915161 ], # C00003
    [  6.60632e-5, 0.000193512 ], # C00004
    [  8.52043e-5, 0.000989574 ], # C00008
    [  0.00564784,        0.02 ], # C00009
    [ 0.000789449,  0.00279558 ], # C00010
    [     1.36e-5,     1.36e-5 ], # C00011
    [  0.00133465,  0.00472625 ], # C00024
    [ 0.000339974, 0.000722382 ], # C00026
    [  1.24527e-7,  1.90545e-6 ], # C00036
    [ 0.000416036,  0.00235402 ], # C00042
#   [    0.000001,         0.1 ], # C00080
    [  7.76413e-5,  0.00331968 ], # C00091
    [ 0.000297244, 0.000512007 ], # C00122
    [  0.00119182,  0.00205293 ], # C00149
    [ 0.000105274,         0.1 ], # C00158
    [  3.67259e-6, 0.000230163 ], # C00417
    [        1e-7, 0.000303877 ], # C00311
    [     0.00005,       0.005 ], # C00390
    [     0.00005,       0.005 ]  # C00399
]))

# Constrain concentration ratios
# Use natural log
ratio_lim = np.log(np.array([
    [ 0.021, 0.19 ], # 0.021 < C00004 / C00003 < 0.19
    [   0.1,   10 ]  # 0.1 < C00399 / C00390 < 10
]))

ratio_mat = np.array([
    [ 0,  0],
    [ 0,  0],
    [-1,  0], # Ratio 0, C00004 as dividend
    [ 1,  0], # Ratio 0, C00003 as divisor
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  0],
    [ 0,  1], # Ratio 1, C00390 as divisor
    [ 0, -1]  # Ratio 1, C00399 as dividend
])

# Define function for random sampling of concentrations
def random_c(c_lim):
    sample = np.array([np.random.random() for n in range(0, c_lim.shape[0])])
    return sample * (c_lim[:,1] - c_lim[:,0]) + c_lim[:,0]

# Define function for checking if set is thermodynamically feasible
def df_ok(c):
    # Calculate delta G prime
    df = -(g + RT * np.sum(np.transpose(S) * c, 1))
    # Check if all driving forces higher than 0
    return sum(df > 0) == df.shape[0]

# Define function for checking if set has acceptable ratios
def ratios_ok(c):
    ratios = np.sum(ratio_mat.T * c, 1).reshape([ratio_lim.shape[0], 1])
    min = np.sum(np.subtract(ratios, ratio_lim) >= 0, 0)[0] == ratios.shape[0]
    max = np.sum(np.subtract(ratios, ratio_lim) <= 0, 0)[1] == ratios.shape[0]
    return min and max

# Define function for checking that sum of concentrations is not too high
def sum_ok(c, max_tot_c = 0.05):
    return np.sum(np.exp(c)) <= max_tot_c

# Define function that checks concentrations are within limits
def limits_ok(c):
    c_l = c.reshape([c.shape[0],1])
    min = np.sum(np.subtract(c_l, c_lim) >= 0, 0)[0] == c.shape[0]
    max = np.sum(np.subtract(c_l, c_lim) <= 0, 0)[1] == c.shape[0]
    return min and max

# Define function for checking feasibility, ratios, sum, and limits in one go
def is_feasible(c):
    return df_ok(c) and ratios_ok(c) and sum_ok(c[2:]) and limits_ok(c)

# Define function for selecting a random direction
def random_direction(c):
    # Create a random vector of the same length as c
    direction = np.array([np.random.random() for n in range(0, c.shape[0])])
    # Subtract 0.5 to introduce negative directions
    direction = direction - 0.5
    # Set fixed concentration direction to zero
    direction[c_lim[:,1] - c_lim[:,0] == 0] = 0
    # Normalize length of direction vector
    normalized_direction = direction / np.linalg.norm(direction)
    return normalized_direction

# Define function to generate one feasible metabolite concentration set
def generate_feasible_c(c_lim):
    c = random_c(c_lim) # Initialize c
    while not is_feasible(c):
        c = random_c(c_lim) # Generate new c until feasible
    return c

# Create concentration datasets where one metabolite straddles the limit
c_lim_extreme_low = np.copy(c_lim)
c_lim_extreme_low[1,1] = c_lim_extreme_low[1,0] # Bind 2nd low
c_extreme_low = generate_feasible_c(c_lim_extreme_low)

c_lim_extreme_high = np.copy(c_lim)
c_lim_extreme_high[5,0] = c_lim_extreme_high[5,1] # Bind 5th high
c_extreme_high = generate_feasible_c(c_lim_extreme_high)

c_lim_extreme_duo_same = np.copy(c_lim)
c_lim_extreme_duo_same[1,0] = c_lim_extreme_duo_same[1,1] # Bind 2nd low
c_lim_extreme_duo_same[5,0] = c_lim_extreme_duo_same[5,1] # Bind 5th high
c_extreme_duo_same = generate_feasible_c(c_lim_extreme_duo_same)

c_lim_extreme_duo_diff = np.copy(c_lim)
c_lim_extreme_duo_diff[1,1] = c_lim_extreme_duo_diff[1,0] # Bind 2nd low
c_lim_extreme_duo_diff[5,0] = c_lim_extreme_duo_diff[5,1] # Bind 5th high
c_extreme_duo_diff = generate_feasible_c(c_lim_extreme_duo_diff)

# Select extreme low or extreme high for testing
c = np.copy(c_extreme_low)
c = np.copy(c_extreme_high)
c = np.copy(c_extreme_duo_same)
c = np.copy(c_extreme_duo_diff)

# Create a random direction
direction = random_direction(c)

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

# Checking that hard limits are fine
negative = []
positive = []
for n in range(0, 1000):
    cn = generate_feasible_c(c_lim)
    dn = random_direction(cn)
    dn = unstick_direction(cn, dn, c_lim)
    theta = calculate_theta_hard_limit(cn, dn, c_lim)
    negative.append(
        not limits_ok(cn + dn * (theta[0] - 1e-7)) and
        not limits_ok(cn + dn * (theta[1] + 1e-7))
    )
    positive.append(
        limits_ok(cn + dn * theta[0]) and limits_ok(cn + dn * theta[1])
    )
# Seems to be ok

# Define function for determining minimum and maximum step length (theta)
def theta_range(c, direction, precision=1e-3):
    # Define function for honing in on a theta limit
    def hone_theta(theta_outer, theta_inner=0):
        if is_feasible(c + theta_outer * direction):
            # If the outer theta is feasible, accept that solution
            theta_inner = theta_outer
        else:
            while abs(theta_outer - theta_inner) > precision:
                # Calculate a theta between outer and inner limits
                theta_cur = (theta_outer + theta_inner) / 2
                if is_feasible(c + theta_cur * direction):
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
def hit_and_run(S, g, c_lim, ratio_lim, ratio_mat, n_samples, precision=1e-3):
    # Generate starting point
    c = generate_feasible_c(c_lim)
    # Set up concentration storage list
    fMCSs = [c]
    # Perform n steps
    for i in range(0, n_samples - 1):
        # Generate random direction
        direction = random_direction(c)
        # Determine minimum and maximum step length
        theta = theta_range(c, direction, precision=precision)
        # Perform a random sampling of the step length
        theta = theta[0] + np.random.random() * (theta[1] - theta[0])
        # Perform step
        c = c + theta * direction
        # Ensure feasibility
        if not is_feasible(c):
            print("Warning: Infeasible point reached.")
            break
        # Store concentration
        fMCSs.append(c)
    # Return list of concentrations
    return fMCSs

# Define function for performing rejection sampling
def rejection_sampling(S, g, c_lim, ratio_lim, ratio_mat, n_samples):
    return [generate_feasible_c(c_lim) for i in range(0, n_samples)]

# Compare speed of both approaches
# import time
#
# start = time.time()
# fMCSs_hitandrun = hit_and_run(S, g, c_lim, ratio_lim, ratio_mat, 10000)
# hitandrun_time = time.time() - start # 7.64 s
#
# start = time.time()
# fMCSs_rejection = rejection_sampling(S, g, c_lim, ratio_lim, ratio_mat, 10000)
# rejection_time = time.time() - start # 33.9 s
#
# start = time.time()
# fMCSs_hitandrun_high_precision = hit_and_run(
#     S, g, c_lim, ratio_lim, ratio_mat, 10000, precision=1e-6
# )
# hitandrun_high_precision_time = time.time() - start # 12.2 s

# Compare results of both approaches
# for i in range(0, 5):
#     for c in hit_and_run(S, g, c_lim, ratio_lim, ratio_mat, 20000):
#         # Print CSV in mM
#         print("hit_and_run," + ",".join([str(np.exp(x)*1000) for x in c]))

for c in hit_and_run(S, g, c_lim, ratio_lim, ratio_mat, 200000):
    # Print CSV in mM
    print("hit_and_run," + ",".join([str(np.exp(x)*1000) for x in c]))

for c in rejection_sampling(S, g, c_lim, ratio_lim, ratio_mat, 200000):
    print("rejection," + ",".join([str(np.exp(x)*1000) for x in c]))
