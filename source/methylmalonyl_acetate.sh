# Perform MDF analysis
ls data/mmcoa_fixed_acetate/* | while read conc_file; do
  # Extract group name (extracellular acetate concentration type)
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
    --outfile results/mmcoa_fixed_acetate/methylmalonyl.${1}.${i}.tab \
    --constraints data/mmcoa_fixed_acetate/methylmalonyl.concentrations.${1}.tab \
    --concs results/methylmalonyl.mdf.${1}.csv \
    --steps 10000000 --max_conc 1.3 -n 1000"
  done | parallel --no-notice --jobs 10 > /dev/null 2>&1 &
}

# Call function to do sampling
run_sampling "ac_0_25x" # o
run_sampling "ac_0_5x" # d
run_sampling "ac_1x" # d
run_sampling "ac_2x" # m
run_sampling "ac_free" # m

# Plot results
./plot_samples.R \
-i results/mmcoa_fixed_acetate \
-S data/methylmalonyl.stoich.tab \
-G data/methylmalonyl.model_drGs.tab \
-c data/mmcoa_fixed_acetate/methylmalonyl.concentrations.ac_free.tab \
-x H2O \
-o results/mmcoa_fixed_acetate_plots

# Plot results with fewer reactions, metabolites, and conditions
./plot_samples.R \
-i results/mmcoa_fixed_acetate \
-S data/methylmalonyl.stoich.tab \
-G data/methylmalonyl.model_drGs.tab \
-c data/mmcoa_fixed_acetate/methylmalonyl.concentrations.ac_free.tab \
-C data/methylmalonyl.plot_config_less.tab \
-o results/mmcoa_fixed_acetate_plots_reduced
