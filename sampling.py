# Import libraries
import numpy as np

# Define function for random sampling of concentrations
def random_c(c_lim):
    sample = np.array([np.random.random() for n in range(0, c_lim.shape[0])])
    return sample * (c_lim[:,1] - c_lim[:,0]) + c_lim[:,0]

# Define function for checking if set is thermodynamically feasible
def df_ok(c, g, RT, S):
    # Calculate delta G prime
    df = -(g + RT * np.sum(np.transpose(S) * c, 1))
    # Check if all driving forces higher than 0
    return sum(df > 0) == df.shape[0]

# Define function for checking if set has acceptable ratios
def ratios_ok(c, ratio_lim, ratio_mat):
    ratios = np.sum(ratio_mat.T * c, 1).reshape([ratio_lim.shape[0], 1])
    min = np.sum(np.subtract(ratios, ratio_lim) > 0, 0)[0] == ratios.shape[0]
    max = np.sum(np.subtract(ratios, ratio_lim) < 0, 0)[1] == ratios.shape[0]
    return min and max

# Define function for checking that sum of concentrations is not too high
def sum_ok(c, max_tot_c):
    return np.sum(np.exp(c)) <= max_tot_c

# Define function that checks concentrations are within limits
def limits_ok(c, c_lim):
    c_l = c.reshape([c.shape[0],1])
    min = np.sum(np.subtract(c_l, c_lim) >= 0, 0)[0] == c.shape[0]
    max = np.sum(np.subtract(c_l, c_lim) <= 0, 0)[1] == c.shape[0]
    return min and max

# Define function for checking feasibility, ratios, sum, and limits in one go
def is_feasible(c, g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim):
    # Run all feasibility tests
    ok = [
        # Driving forces positive
        df_ok(c, g, RT, S),
        # Concentration ratios within bounds
        ratios_ok(c, ratio_lim, ratio_mat),
        # Concentration sum below max
        sum_ok(c, max_tot_c),
        # Concentrations within bounds
        limits_ok(c, c_lim)
    ]
    # Check that all feasibility tests are True
    return sum(ok) == len(ok)

# Define function for selecting a random direction
def random_direction(c, c_lim):
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
def generate_feasible_c(g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim):
    c = random_c(c_lim) # Initialize c
    while not is_feasible(c, g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim):
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
    c_lim, direction, precision=1e-3
    ):
    # Define function for honing in on a theta limit
    def hone_theta(theta_outer, theta_inner=0):
        if is_feasible(
            c + theta_outer * direction, g, RT, S, ratio_lim, ratio_mat,
            max_tot_c, c_lim
            ):
            # If the outer theta is feasible, accept that solution
            theta_inner = theta_outer
        else:
            while abs(theta_outer - theta_inner) > precision:
                # Calculate a theta between outer and inner limits
                theta_cur = (theta_outer + theta_inner) / 2
                if is_feasible(
                    c + theta_cur * direction, g, RT, S, ratio_lim, ratio_mat,
                    max_tot_c, c_lim
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
    c_lim, direction, n_samples, precision=1e-3
    ):
    # Set up concentration storage list
    fMCSs = [c]
    # Perform n steps
    for i in range(0, n_samples - 1):
        # Generate random direction
        direction = random_direction(c, c_lim)
        # Determine minimum and maximum step length
        theta = theta_range(
            c, g, RT, S, ratio_lim, ratio_mat, max_tot_c,
            c_lim, direction, precision=1e-3
        )
        # Perform a random sampling of the step length
        theta = theta[0] + np.random.random() * (theta[1] - theta[0])
        # Perform step
        c = c + theta * direction
        # Ensure feasibility
        if not is_feasible(c, g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim):
            print("Warning: Infeasible point reached.")
            break
        # Store concentration
        fMCSs.append(c)
    # Return list of concentrations
    return fMCSs
