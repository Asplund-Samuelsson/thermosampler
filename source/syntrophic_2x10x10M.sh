#!/bin/bash
# Create output directory
mkdir -p results/syntrophic_2x10x10M

# Sample 10 replicates of 10 million steps in parallel at feasible and optimum delta G's
for i in {1..10}; do
  echo "./sampling.py --reactions data/syntrophic.model.tab \
  --std_drG data/syntrophic.model_drGs.tab \
  --ratios data/syntrophic.ratio_range.tab \
  --outfile results/syntrophic_2x10x10M/syntrophic.sampling.${i}.tab \
  --constraints data/syntrophic.concentrations.tab \
  --concs results/syntrophic_mdf.csv --steps 10000000 --max_conc 1.5 \
  --conc_sums data/syntrophic.conc_sums.tab -n 1000 --quiet"

  echo "./sampling.py --reactions data/syntrophic.model.tab \
  --std_drG data/syntrophic.model_drGs.tab \
  --ratios data/syntrophic.ratio_range.tab \
  --outfile results/syntrophic_2x10x10M/syntrophic.mdf_sampling.${i}.tab \
  --constraints data/syntrophic.concentrations.tab \
  --concs results/syntrophic_mdf.csv --steps 10000000 --max_conc 1.5 \
  --mdf 0.447935 --conc_sums data/syntrophic.conc_sums.tab -n 1000 --quiet"
done | parallel --no-notice --jobs -2 --bar
