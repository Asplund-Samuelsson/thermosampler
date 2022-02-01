#!/bin/bash
# Create output directory
mkdir -p results/methylmalonyl_2x10x1M

# Sample 10 replicates of 1 million steps in parallel at feasible and optimum delta G's
for i in {1..10}; do
  echo "./sampling.py --reactions data/methylmalonyl.model.tab \
  --std_drG data/methylmalonyl.model_drGs.tab \
  --ratios data/methylmalonyl.ratio_range.tab \
  --outfile results/methylmalonyl_2x10x1M/methylmalonyl.sampling.${i}.tab \
  --constraints data/methylmalonyl.concentrations.tab \
  --concs results/methylmalonyl_mdf.csv --steps 1000000 --max_conc 1.1 --quiet"

  echo "./sampling.py --reactions data/methylmalonyl.model.tab \
  --std_drG data/methylmalonyl.model_drGs.tab \
  --ratios data/methylmalonyl.ratio_range.tab \
  --outfile results/methylmalonyl_2x10x1M/methylmalonyl.mdf_sampling.${i}.tab \
  --constraints data/methylmalonyl.concentrations.tab \
  --concs results/methylmalonyl_mdf.csv --steps 1000000 --max_conc 1.1 \
  --mdf 1.480283 --quiet"
done | parallel --no-notice --jobs -2 --bar
