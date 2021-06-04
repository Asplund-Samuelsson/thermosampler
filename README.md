# Thermosampler

Linear programming and Markov Chain Monte Carlo sampling tools for thermodynamic analysis of metabolism.

### MDF and NEM analysis

`mdf.py` performs Max-min Driving Force (MDF; [Noor _et al._, 2014](http://doi.org/10.1371/journal.pcbi.1003483)) and Network-Embedded
MDF (NEM) analysis.

##### _Example: Perform NEM analysis on lysine biosynthesis in_ E. coli _and_ Synechocystis

`./mdf.py --min_conc 0.0000001 --max_conc 0.1 --constraints examples/E_coli.concentrations.tab --ratios examples/E_coli.Lys_opt_ratios.tab --pathway examples/E_coli.Lys_pathway.txt examples/E_coli.model.tab examples/E_coli.model_drgs.tab results/E_coli_Lys_nem.csv`

`./mdf.py --min_conc 0.0000001 --max_conc 0.1 --constraints examples/Synechocystis.concentrations.tab --ratios examples/Synechocystis.Lys_opt_ratios.tab --pathway examples/Synechocystis.Lys_pathway.txt examples/Synechocystis.model.tab examples/Synechocystis.model_drgs.tab results/Synechocystis_Lys_nem.csv`

##### _Example: Perform MDF analysis on TCA cycle_

`./mdf.py --min_conc 0.0000001 --max_conc 0.1 --constraints examples/tca.concentrations.tab --ratios examples/tca.ratios.tab examples/tca.model.tab examples/tca.model_drgs.tab results/tca_mdf.csv`
