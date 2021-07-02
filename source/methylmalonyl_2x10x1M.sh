# Sample 10 replicates of 1 million steps in parallel at feasible delta G
for i in {1..10}; do
  echo "./sampling.py --reactions data/methylmalonyl.model.tab \
  --std_drG data/methylmalonyl.model_drGs.tab \
  --ratios data/methylmalonyl.ratio_range.tab \
  --outfile results/methylmalonyl_2x10x1M/methylmalonyl.sampling.${i}.tab \
  --constraints data/methylmalonyl.concentrations.tab \
  --concs results/methylmalonyl_mdf.csv --steps 1000000 --max_conc 1.1"
done | parallel --no-notice --jobs 10

# Sample 10 replicates of 1 million steps in parallel at optimum delta G
for i in {1..10}; do
  echo "./sampling.py --reactions data/methylmalonyl.model.tab \
  --std_drG data/methylmalonyl.model_drGs.tab \
  --ratios data/methylmalonyl.ratio_range.tab \
  --outfile results/methylmalonyl_2x10x1M/methylmalonyl.mdf_sampling.${i}.tab \
  --constraints data/methylmalonyl.concentrations.tab \
  --concs results/methylmalonyl_mdf.csv --steps 1000000 --max_conc 1.1 \
  --mdf 1.480283"
done | parallel --no-notice --jobs 10
