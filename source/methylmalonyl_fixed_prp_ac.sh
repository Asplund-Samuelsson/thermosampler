#!/bin/bash
# Create output directory
mkdir -p results/mmcoa_fixed_prp_ac

# Perform MDF analysis
ls data/mmcoa_fixed_prp_ac/* | while read conc_file; do
  # Extract group name (extracellular acetate and propionate concentration type)
  Group=`echo $conc_file | cut -f 3- -d . | sed -e 's/.tab//'`
  # Do MDF
  ./mdf.py --min_conc 1e-6 --max_conc 1e-2 \
  --constraints $conc_file \
  --ratios data/methylmalonyl.ratios.tab \
  data/methylmalonyl.model.tab \
  data/methylmalonyl.model_drGs.tab \
  results/methylmalonyl.mdf.$Group.csv
done

# Sample 10 replicates of 10 million steps in parallel at feasible delta G
run_sampling () {
  for i in {1..10}; do
    echo "./sampling.py --reactions data/methylmalonyl.model.tab \
    --std_drG data/methylmalonyl.model_drGs.tab \
    --ratios data/methylmalonyl.ratio_range.tab \
    --outfile results/mmcoa_fixed_prp_ac/methylmalonyl.${1}.${i}.tab \
    --constraints data/mmcoa_fixed_prp_ac/methylmalonyl.concentrations.${1}.tab \
    --concs results/methylmalonyl.mdf.${1}.csv \
    --steps 1000000 --max_conc 1.3 -n 100"
  done | parallel --no-notice --jobs -2 > /dev/null 2>&1 &
}

# Call function to do sampling
run_sampling "free"
run_sampling "prpLo"
run_sampling "acHi"
run_sampling "prpLo_acHi"
run_sampling "prpHi_acLo"

# Plot results
./plot_samples.R \
-i results/mmcoa_fixed_prp_ac \
-S data/methylmalonyl.stoich.tab \
-G data/methylmalonyl.model_drGs.tab \
-c data/mmcoa_fixed_prp_ac/methylmalonyl.concentrations.free.tab \
-C data/methylmalonyl.plot_config.fixed_prp_ac.tab \
-x H2O \
-o results/mmcoa_fixed_prp_ac_plots

# Plot results with fewer reactions, metabolites, and conditions
./plot_samples.R \
-i results/mmcoa_fixed_prp_ac \
-S data/methylmalonyl.stoich.tab \
-G data/methylmalonyl.model_drGs.tab \
-c data/mmcoa_fixed_prp_ac/methylmalonyl.concentrations.free.tab \
-C data/methylmalonyl.plot_config.fixed_prp_ac_less.tab \
-o results/mmcoa_fixed_prp_ac_plots_reduced
