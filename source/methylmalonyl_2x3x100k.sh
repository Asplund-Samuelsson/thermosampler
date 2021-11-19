#!/bin/bash
# Feasible, Replicate 1
./sampling.py --reactions data/methylmalonyl.model.tab \
--std_drG data/methylmalonyl.model_drGs.tab \
--ratios data/methylmalonyl.ratio_range.tab \
--outfile results/methylmalonyl.sampling.tab \
--constraints data/methylmalonyl.concentrations.tab \
--concs results/methylmalonyl_mdf.csv --steps 100000 --max_conc 1.1

# Optimum, Replicate 1
./sampling.py --reactions data/methylmalonyl.model.tab \
--std_drG data/methylmalonyl.model_drGs.tab \
--ratios data/methylmalonyl.ratio_range.tab \
--outfile results/methylmalonyl.mdf_sampling.tab \
--constraints data/methylmalonyl.concentrations.tab \
--concs results/methylmalonyl_mdf.csv --steps 100000 \
--max_conc 1.1 --mdf 1.480283

# Feasible, Replicate 2
./sampling.py --reactions data/methylmalonyl.model.tab \
--std_drG data/methylmalonyl.model_drGs.tab \
--ratios data/methylmalonyl.ratio_range.tab \
--outfile results/methylmalonyl.sampling_rep2.tab \
--constraints data/methylmalonyl.concentrations.tab \
--concs results/methylmalonyl_mdf.csv --steps 100000 --max_conc 1.1

# Optimum, Replicate 2
./sampling.py --reactions data/methylmalonyl.model.tab \
--std_drG data/methylmalonyl.model_drGs.tab \
--ratios data/methylmalonyl.ratio_range.tab \
--outfile results/methylmalonyl.mdf_sampling_rep2.tab \
--constraints data/methylmalonyl.concentrations.tab \
--concs results/methylmalonyl_mdf.csv --steps 100000 \
--max_conc 1.1 --mdf 1.480283

# Feasible, Replicate 3
./sampling.py --reactions data/methylmalonyl.model.tab \
--std_drG data/methylmalonyl.model_drGs.tab \
--ratios data/methylmalonyl.ratio_range.tab \
--outfile results/methylmalonyl.sampling_rep3.tab \
--constraints data/methylmalonyl.concentrations.tab \
--concs results/methylmalonyl_mdf.csv --steps 100000 --max_conc 1.1

# Optimum, Replicate 3
./sampling.py --reactions data/methylmalonyl.model.tab \
--std_drG data/methylmalonyl.model_drGs.tab \
--ratios data/methylmalonyl.ratio_range.tab \
--outfile results/methylmalonyl.mdf_sampling_rep3.tab \
--constraints data/methylmalonyl.concentrations.tab \
--concs results/methylmalonyl_mdf.csv --steps 100000 \
--max_conc 1.1 --mdf 1.480283
