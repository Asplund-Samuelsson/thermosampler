# Sample 10 replicates of 10 million steps in parallel at feasible delta G
for i in {1..10}; do
  echo "./sampling.py --reactions data/syntrophic.model.tab \
  --std_drG data/syntrophic.model_drGs.tab \
  --ratios data/syntrophic.ratio_range.tab \
  --outfile results/syntrophic_2x10x10M/syntrophic.high_CH4.${i}.tab \
  --constraints data/syntrophic.concentrations.high_CH4.tab \
  --concs results/syntrophic_mdf.high_CH4.csv --steps 10000000 --max_conc 1.5 \
  --conc_sums data/syntrophic.conc_sums.tab -n 1000"
done | parallel --no-notice --jobs 10 > /dev/null 2>&1 &

# Sample 10 replicates of 10 million steps in parallel at optimum delta G
for i in {1..10}; do
  echo "./sampling.py --reactions data/syntrophic.model.tab \
  --std_drG data/syntrophic.model_drGs.tab \
  --ratios data/syntrophic.ratio_range.tab \
  --outfile results/syntrophic_2x10x10M/syntrophic.mdf_high_CH4.${i}.tab \
  --constraints data/syntrophic.concentrations.high_CH4.tab \
  --concs results/syntrophic_mdf.high_CH4.csv --steps 10000000 --max_conc 1.5 \
  --mdf 0.447935 --conc_sums data/syntrophic.conc_sums.tab -n 1000"
done | parallel --no-notice --jobs 10 > /dev/null 2>&1 &

# Plot results including MDF
./plot_samples.R \
-i results/syntrophic_2x10x10M \
-S data/syntrophic.stoich.tab \
-G data/syntrophic.model_drGs.tab \
-c data/syntrophic.concentrations.tab \
-x H2O \
-o results/syntrophic_CH4_mdf

# Plot without MDF data for more clarity
mkdir results/syntrophic_CH4_comparison

ls results/syntrophic_2x10x10M/* | grep -v "mdf" | while read Infile; do
  Outfile=`
    realpath $Infile | sed -e 's/2x10x10M/CH4_comparison/' \
    -e 's/sampling/wide_CH4/'
  `
  ln -s `realpath $Infile` $Outfile
done

./plot_samples.R \
-i results/syntrophic_CH4_comparison \
-S data/syntrophic.stoich.tab \
-G data/syntrophic.model_drGs.tab \
-c data/syntrophic.concentrations.tab \
-x H2O \
-o results/syntrophic_CH4
