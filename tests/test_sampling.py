#!/usr/bin/env python3

# Add repository root to the path
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

# Import the script to be tested
from sampling import *

# Ensure testing is performed on the same random seed every time
np.random.seed(42)

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
    [-1,  0], # Ratio 0, C00003 as divisor
    [ 1,  0], # Ratio 0, C00004 as dividend
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
    [ 0, -1], # Ratio 1, C00390 as divisor
    [ 0,  1]  # Ratio 1, C00399 as dividend
])

# Select maximum sum of concentrations
max_tot_c = 0.05 + 1 # Add 1 M for water

# MDF concentration set
c_mdf = np.log(np.array([
    1.0000000000, # C00001
    0.0007687213, # C00002
    0.0018026464, # C00003
    0.0001135667, # C00004
    0.0004237077, # C00008
    0.0118182519, # C00009
    0.0015460872, # C00010
    0.0000136000, # C00011
    0.0025373290, # C00024
    0.0004954426, # C00026
    0.0000001614, # C00036
    0.0006820748, # C00042
    0.0014036845, # C00091
    0.0005120070, # C00122
    0.0011918200, # C00149
    0.0078406055, # C00158
    0.0000243520, # C00417
    0.0000018146, # C00311
    0.0005000000, # C00390
    0.0005000000  # C00399
]))

# Define tests
def test_random_c(c_lim=c_lim):
    # Sample 1000 times and make sure the data is new and not out of bounds
    prev = random_c(c_lim)
    for i in range(0,1000):
        # Pick random concentrations
        random_concentrations = random_c(c_lim)
        # Check that they are not out of bounds
        assert not False in (random_concentrations >= c_lim[:,0])
        assert not False in (random_concentrations <= c_lim[:,1])
        # Check that they are not stuck
        assert not np.array_equal(random_concentrations, c_lim[:,0])
        assert not np.array_equal(random_concentrations, c_lim[:,1])
        # Check that they change between rounds
        assert not np.array_equal(prev, random_concentrations)
        prev = random_concentrations.copy()

def test_df_ok(RT=RT):
    S_small = np.array(
        [[-1,  0],
         [ 1, -1],
         [ 0,  1]]
    )
    g_small = np.array([-10, 10])
    c_small = np.log(np.array([0.002, 0.1, 0.0015]))
    c_small_fail = np.log(np.array([0.002, 0.1005, 0.002]))
    assert df_ok(c_small, g_small, RT, S_small)
    assert not df_ok(c_small, np.array(list(reversed(g_small))), RT, S_small)
    assert not df_ok(c_small_fail, g_small, RT, S_small)

def test_ratios_ok(ratio_lim=ratio_lim, ratio_mat=ratio_mat):
    # Check ratios in middle
    c_rat = np.array([
        0, 0,
        np.log(0.01), np.log(0.001),
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        np.log(0.05), np.log(0.01)
    ])
    assert ratios_ok(c_rat, ratio_lim, ratio_mat)
    # Check ratios near limits
    c_rat = np.array([
        0, 0,
        np.log(0.01), np.log(0.00189999),
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        np.log(0.199999), np.log(0.02)
    ])
    assert ratios_ok(c_rat, ratio_lim, ratio_mat)
    # Check ratios at limits
    c_rat = np.array([
        0, 0,
        np.log(0.01), np.log(0.0019),
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        np.log(0.2), np.log(0.02)
    ])
    assert not ratios_ok(c_rat, ratio_lim, ratio_mat)
    # Check one ratio beyond limits
    c_rat = np.array([
        0, 0,
        np.log(0.002), np.log(0.01),
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        np.log(0.01), np.log(0.05)
    ])
    assert not ratios_ok(c_rat, ratio_lim, ratio_mat)
    # Check other ratio beyond limits
    c_rat = np.array([
        0, 0,
        np.log(0.001), np.log(0.01),
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        np.log(0.1), np.log(0.005)
    ])
    assert not ratios_ok(c_rat, ratio_lim, ratio_mat)
    # Check both ratios beyond limits
    c_rat = np.array([
        0, 0,
        np.log(0.0002), np.log(0.01),
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        np.log(0.01), np.log(0.2)
    ])
    assert not ratios_ok(c_rat, ratio_lim, ratio_mat)

def test_sum_ok():
    c_success = np.log(np.array([0.01,0.01,0.01]))
    c_fail = np.log(np.array([0.0001,0.03,0.01,0.1]))
    assert sum_ok(c_success, 0.04)
    assert sum_ok(c_success, 0.03000001)
    assert not sum_ok(c_success, 0.01)
    assert not sum_ok(c_fail, 0.04)
    assert sum_ok(np.array([0,-3]), 0.5 + 1)

def test_limits_ok():
    c_lim_test = np.log(np.array(
        [[0.01,  0.1  ],
         [0.001, 0.01 ],
         [0.002, 0.002]]
    ))
    c_pass_1 = np.log(np.array([0.05, 0.005, 0.002]))
    c_pass_2 = np.log(np.array([0.1, 0.001, 0.002]))
    c_fail_1 = np.log(np.array([0.009, 0.005, 0.002]))
    c_fail_2 = np.log(np.array([0.05, 0.02, 0.002]))
    c_fail_3 = np.log(np.array([0.05, 0.005, 0.0021]))
    assert limits_ok(c_pass_1, c_lim_test)
    assert limits_ok(c_pass_2, c_lim_test)
    assert not limits_ok(c_fail_1, c_lim_test)
    assert not limits_ok(c_fail_2, c_lim_test)
    assert not limits_ok(c_fail_3, c_lim_test)

def test_is_feasible(RT=RT):
    S_small = np.array(
        [[-1,  1],
         [-1,  0],
         [ 1, -1],
         [ 0,  1]]
    )
    g_small = np.array([-10, 10])
    c_small = np.log(np.array([1.0, 0.002, 0.1, 0.0015]))
    # Violates dG
    c_small_fail_1 = np.log(np.array([1.0, 0.002, 0.1005, 0.002]))
    # Violates ratio
    c_small_fail_2 = np.log(np.array([1.0, 0.0031, 0.1, 0.0015]))
    # Violates total conc
    c_small_fail_3 = np.log(np.array([1.0, 0.02, 0.2, 0.002]))
    # Violates limits
    c_small_fail_4 = np.log(np.array([1.0, 0.002, 0.1, 0.0001]))
    # Violates everything
    c_small_fail_5 = np.log(np.array([1.0, 0.00001, 0.5, 0.5]))
    ratio_lim_small = np.log(np.array([[0.01, 0.03]]))
    ratio_mat_small = np.array([[0,1,-1,0]]).transpose()
    max_tot_c_test = 0.11 + 1
    c_lim_test = np.log(np.array([
        [1.0, 1.0],
        [0.001, 0.01],
        [0.05, 0.2],
        [0.001, 0.03]
    ]))
    assert is_feasible(
        c_small, g_small, RT, S_small, ratio_lim_small, ratio_mat_small,
        max_tot_c_test, c_lim_test
    )
    assert not is_feasible(
        c_small_fail_1, g_small, RT, S_small, ratio_lim_small, ratio_mat_small,
        max_tot_c_test, c_lim_test
    )
    assert not is_feasible(
        c_small_fail_2, g_small, RT, S_small, ratio_lim_small, ratio_mat_small,
        max_tot_c_test, c_lim_test
    )
    assert not is_feasible(
        c_small_fail_3, g_small, RT, S_small, ratio_lim_small, ratio_mat_small,
        max_tot_c_test, c_lim_test
    )
    assert not is_feasible(
        c_small_fail_4, g_small, RT, S_small, ratio_lim_small, ratio_mat_small,
        max_tot_c_test, c_lim_test
    )
    assert not is_feasible(
        c_small_fail_5, g_small, RT, S_small, ratio_lim_small, ratio_mat_small,
        max_tot_c_test, c_lim_test
    )

def test_random_direction():
    # Sample 1000 times
    prev = random_direction(c_lim)
    for i in range(0,1000):
        # Pick random direction
        random_dir = random_direction(c_lim)
        # Check direction is zero for fixed concentrations
        assert not np.count_nonzero(random_dir[c_lim[:,0] == c_lim[:,1]])
        # Check that they are not stuck
        assert not np.array_equal(random_dir, c_lim)
        # Check that they change between rounds
        assert not np.array_equal(prev, random_dir)
        prev = random_dir.copy()

def test_generate_feasible_c():
    prev = generate_feasible_c(g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim)
    for i in range(0,1000):
        # Generate feasible concentrations
        c = generate_feasible_c(
            g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim
        )
        # Make sure result is feasible
        assert is_feasible(c, g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim)
        # Check that they change between rounds
        assert not np.array_equal(prev, c)
        # Make copy for comparison in next round
        prev = c.copy()

def test_unstick_direction():
    # Checking that hard limits are fine
    negative = []
    positive = []
    for n in range(0, 1000):
        cn = generate_feasible_c(
            g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim
        )
        dn = random_direction(c_lim)
        dn = unstick_direction(cn, dn, c_lim)
        theta = calculate_theta_hard_limit(cn, dn, c_lim)
        negative.append(
            not limits_ok(cn + dn * (theta[0] - 1e-7), c_lim) and
            not limits_ok(cn + dn * (theta[1] + 1e-7), c_lim)
        )
        positive.append(
            limits_ok(cn + dn * theta[0], c_lim) and \
            limits_ok(cn + dn * theta[1], c_lim)
        )
    # Assert that negative and positive are fine
    assert sum(negative) == len(negative)
    assert sum(positive) == len(positive)

    for i in range(0, 100):

        # Create concentration datasets where one metabolite straddles the limit
        c_lim_extreme_low = np.copy(c_lim)
        c_lim_extreme_low[1,1] = c_lim_extreme_low[1,0] # Bind 2nd low
        c_extreme_low = generate_feasible_c(
            g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim_extreme_low
        )

        c_lim_extreme_high = np.copy(c_lim)
        c_lim_extreme_high[5,0] = c_lim_extreme_high[5,1] # Bind 6th high
        c_extreme_high = generate_feasible_c(
            g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim_extreme_high
        )

        c_lim_extreme_duo_same = np.copy(c_lim)
        c_lim_extreme_duo_same[2,0] = c_lim_extreme_duo_same[2,1] # Bind 3rd high
        c_lim_extreme_duo_same[5,0] = c_lim_extreme_duo_same[5,1] # Bind 6th high
        c_extreme_duo_same = generate_feasible_c(
            g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim_extreme_duo_same
        )

        c_lim_extreme_duo_diff = np.copy(c_lim)
        c_lim_extreme_duo_diff[1,1] = c_lim_extreme_duo_diff[1,0] # Bind 2nd low
        c_lim_extreme_duo_diff[5,0] = c_lim_extreme_duo_diff[5,1] # Bind 6th high
        c_extreme_duo_diff = generate_feasible_c(
            g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim_extreme_duo_diff
        )

        # Create a random direction
        direction = random_direction(c_lim)

        # Ensure a small step is possible in one direction
        def small_step(c):
            unstuck_direction = unstick_direction(c, direction, c_lim)
            dir_pos = is_feasible(
                c + 1e-9*unstuck_direction,
                g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim
            )
            dir_neg = is_feasible(
                c - 1e-9*unstuck_direction,
                g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim
            )
            return dir_pos or dir_neg

        assert small_step(c_extreme_low)
        assert small_step(c_extreme_high)
        assert small_step(c_extreme_duo_same)
        assert small_step(c_extreme_duo_diff)

def test_calculate_theta_hard_limit():
    for i in range(0, 1000):
        c = generate_feasible_c(
            g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim
        )
        direction = random_direction(c_lim)
        hard_theta = calculate_theta_hard_limit(c, direction, c_lim)
        # At limits should be OK
        assert limits_ok(c + hard_theta[0]*direction, c_lim)
        assert limits_ok(c + hard_theta[1]*direction, c_lim)
        # Within limits should be OK
        assert limits_ok(
            c + (hard_theta[0] + 1e-6*np.abs(hard_theta[0]))*direction, c_lim
        )
        assert limits_ok(
            c + (hard_theta[1] - 1e-6*np.abs(hard_theta[1]))*direction, c_lim
        )
        # Beyond limits should NOT be OK
        assert not limits_ok(
            c + (hard_theta[0] - 1e-6*np.abs(hard_theta[0]))*direction, c_lim
        )
        assert not limits_ok(
            c + (hard_theta[1] + 1e-6*np.abs(hard_theta[1]))*direction, c_lim
        )

def test_theta_range():
    for i in range(0, 1000):
        c = generate_feasible_c(
            g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim
        )
        direction = random_direction(c_lim)
        t = theta_range(
            c, g, RT, S, ratio_lim, ratio_mat, max_tot_c,
            c_lim, direction, precision=1e-3
        )
        for x in np.linspace(t[0], t[1], 50):
            assert is_feasible(
                c + x*direction,
                g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim
            )

def test_hit_and_run():
    # Check that MDF concentrations are feasible
    assert is_feasible(c_mdf, g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim)
    # Get feasible metabolite concentration sets starting from MDF
    fMCSs = hit_and_run(
        c_mdf, g, RT, S, ratio_lim, ratio_mat, max_tot_c,
        c_lim, n_samples=50000, precision=1e-3
    )
    # Check that each fMCS is feasible
    for i in range(len(fMCSs)):
        assert is_feasible(
            fMCSs[i],
            g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim
        )
        # Check that the fMCSs are different
        if i > 0:
            assert not np.array_equal(fMCSs[i-1], fMCSs[i])
    # Generate fMCSs from several starting points
    for i in range(0, 1000):
        c = generate_feasible_c(
            g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim
        )
        fMCSs = hit_and_run(
            c, g, RT, S, ratio_lim, ratio_mat, max_tot_c,
            c_lim, n_samples=50, precision=1e-3
        )
        # Check that all fMCSs are feasible
        for i in range(len(fMCSs)):
            assert is_feasible(
                fMCSs[i],
                g, RT, S, ratio_lim, ratio_mat, max_tot_c, c_lim
            )
            # Check that the fMCSs are different
            if i > 0:
                assert not np.array_equal(fMCSs[i-1], fMCSs[i])
