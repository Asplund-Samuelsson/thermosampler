![alt text](thermosampler.png "Thermosampler")

# Thermosampler

Linear programming and Markov Chain Monte Carlo sampling tools for thermodynamic analysis of metabolism.

### Contents

1. [Background](#background)
2. [System requirements](#requirements)
3. [MDF and NEM analysis](#mdf)
4. [Hit-and-run analysis](#hit)
5. [Case study: Syntrophic communities](#syntrophism)
6. [Author](#author)

<a name="background"></a>
## Background

### Motivation

Tools such as NET analysis ([Kümmel _et al._ 2006](https://doi.org/10.1038/msb4100074), [Zamboni _et al._ 2008](https://doi.org/10.1186/1471-2105-9-199)) allow calculation of the minimum and maximum feasible metabolite concentrations in a metabolic network, while ensuring that all reactions have a negative Gibbs free energy change, or, equivalently, a positive thermodynamic driving force. Similarly, Max-min Driving Force analysis (MDF; [Noor _et al._, 2014](http://doi.org/10.1371/journal.pcbi.1003483)) is another linear programming approach to finding an optimal set of metabolite concentrations in a metabolic network. However, the aforementioned approaches only identify extreme concentration values, which is why this `thermosampler` framework was developed to explore the full "solution space" of feasible metabolite concentrations (**[Fig 1](#fig1)**). A random walk, or hit-and-run algorithm, allows assessment of metabolite concentration distributions through exploration of all combinations of concentrations that are thermodynamically feasible.

<a name="fig1"></a>

| ![alt text](examples/images/thermosampler_explanation.png "Explanation of the thermosampler hit-and-run algorithm") |
| --- |
| **Fig 1. Thermosampler framework overview and example.** One feasible combination of concentrations within the solution space serves as starting point ("Start") for a random walk. First, a random direction is picked (1). Stepping outside and back into feasible space over and over (2) determines the longest possible step size in the positive and negative directions (3). Then a random step length is performed (4). A new point in the solution space is thus reached and the process is repeated many times over. The right panel depicts 200 steps performed by the `sampling.py` script. Color indicates step number (fMCS; feasible metabolite concentration set). Note that the random walk occurs in logarithmic space, which makes it appear as if higher concentrations are more sparsely sampled. |

### Toy model

As an example, a clearly defined solution space was devised using a "toy" model. The toy model has two metabolites, A and B, and one reaction where A is transformed into B. The standard Gibbs free energy change was selected to be 5.705 kJ/mol, which corresponds to the constraint that [A] must be at least 10 times higher than [B]. Furthermore, [A] and [B] must be within the concentration range 0.01 to 1, and the sum of [A] and [B] must be 1 or lower. These three constraints define a solution space on the AB-plane that has the shape of a triangle (**[Fig 1](#fig1)**). A real metabolic network has dozens of reactions and metabolites, and is therefore many magnitudes more complex and difficult to visualize.

### Example sampling

A random walk applying the `thermosampler` algorithm was performed within the triangular AB solution space using the `sampling.py` script:
```
./sampling.py --reactions examples/toy.model.tab \
--std_drG examples/toy.model_drgs.tab \
--constraints examples/toy.concentrations.tab \
--max_conc 1 --steps 200 \
--outfile results/toy.sampling.tab
```

The results from the random walk were then plotted to yield the right panel in **[Fig 1](#fig1)**:
```
examples/plot_toy_example.R
results/toy.sampling.pdf # Output PDF
```

<a name="feasibility"></a>
### Feasibility checks

In the `thermosampler` algorithm, the _feasibility_ concept and how to check it is very important, as it ensures that the random walk occurs within the solution space, and thereby defines the solution space itself. _Feasibility_ consists of adherence of the metabolite concentration vector to several constraints:
1. **Thermodynamic driving forces** of all reactions must be positive, based on the user-supplied model reactions and directions, as well as the user-supplied standard reaction Gibbs free energy changes.
2. **Metabolite concentration ratios** must be within the user-supplied ranges.
3. **Total metabolite concentration** sum must be below the user-defined value.
4. **Metabolite concentrations** must be within the user-supplied ranges.
5. **Sums of groups of metabolite concentrations** must be below the user-supplied values.

Several functions exist within `sampling.py` to perform the feasibility checks, and are called by the "master" feasibility function named `is_feasible`. It should be relatively simple to add additional feasibility checks to constrain the solution space further.

### Algorithm

The `thermosampler` algorithm starts at a feasible combination of concentrations within the solution space, _i.e._ one "feasible metabolite concentration set" (fMCS). This set is obtained through inefficient rejection sampling or from MDF analysis. More fMCSs are obtained along the random walk.

The `thermosampler` algorithm works as follows, with concentrations in logarithmic space:

1. **Create a step direction vector** of the same length as the concentrations vector, but with random values.
2. **Hone in on the edge of the solution space** by determining the maximum and minimum allowable factor _theta_ to apply to the step direction vector before adding it to the current concentrations vector. First the maximum step size without violating any of the initial concentration range bounds is added and [_feasibility_](#feasibility) is checked. If feasible, the edge has been found. If not, it is necessary to hone in on the edge of solution space by iteratively increasing the inner step size (starts at zero) and decreasing the outer step size (starts at the previously mentioned maximum), until the difference is smaller than a user-defined, small value.
3. Once the difference between inner and outer step size is small enough, the inner step size is taken as "touching" the edge of the solution space, and thereby constituting the maximum feasible step size. The minimum feasible step size (can be negative) is gained from performing step 2 in the opposite direction. Thereby **the minimum and maximum step size is determined**.
4. **Perform one step in the random walk** by sampling one step size factor _theta_ from the previously calculated range, multiplying the step direction vector by that factor, and adding the product to the current metabolite concentrations vector. The feasibility of the new vector, _i.e._ the new metabolite concentration set, is checked. The algorithm stops and prints an error if an infeasible point is reached, which could happen for example if the solution space would be noncontinuous or concave. The error is however more likely to occur from infeasible model parameters (the input files and variables).
5. **Repeat from step 1** until the desired number of steps through the solution space have been performed.

Practical examples of using the `thermosampler` algorithm with `sampling.py` are [shown below](#hit).

<a name="requirements"></a>
## System requirements

- Operating system with Unix-like terminal (Tested on Ubuntu 18.04.5 LTS, 20.04.2 LTS, and MacOS 12.1)
- bash (Tested with 3.2.57, 4.4.20 and 5.0.17)
- Python ≥ 3.7 (Tested with 3.7.6 and 3.10.0)
- R ≥ 4.1.1 (Tested with 4.1.1)
- GNU parallel (Tested with 20161222 and 20220122)
- Python libraries: numpy, pandas, scipy, tqdm, joblib
- R libraries: doMC, foreach, ggridges, optparse, scales, tidyverse

### Setting up the virtual environment

It is recommended to run the Python scripts within a virtual environment. The virtual environment can be created using the following command:

```
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

On Windows systems, replace "python3" with "py", and activate the virtual environment using the command ```venv\Scripts\activate.bat```. Please note that this software has not been tested on Windows systems.

For improved performance om Mac systems with M1 processors, numpy can be built with support for the Accelerate framework. This can lead to roughly a doubling of performance in hit-and-run sampling. For more information, see this [Stack Overflow question](https://stackoverflow.com/questions/69848969/how-to-build-numpy-from-source-linked-to-apple-accelerate-framework). 

<a name="mdf"></a>
## MDF and NEM analysis

`mdf.py` performs Max-min Driving Force (MDF; [Noor _et al._, 2014](http://doi.org/10.1371/journal.pcbi.1003483)) and Network-Embedded
MDF (NEM; [Asplund-Samuelsson _et al._, 2018](https://doi.org/10.1016/j.ymben.2017.12.011)) analysis. Run `./mdf.py -h` to list all options of the script.

#### Example: NEM analysis on lysine biosynthesis in _E. coli_ and _Synechocystis_

```
./mdf.py --min_conc 0.0000001 --max_conc 0.1 \
--constraints examples/E_coli.concentrations.tab \
--ratios examples/E_coli.Lys_opt_ratios.tab \
--pathway examples/E_coli.Lys_pathway.txt \
examples/E_coli.model.tab \
examples/E_coli.model_drgs.tab \
results/E_coli_Lys_nem.csv
```

```
./mdf.py --min_conc 0.0000001 --max_conc 0.1 \
--constraints examples/Synechocystis.concentrations.tab \
--ratios examples/Synechocystis.Lys_opt_ratios.tab \
--pathway examples/Synechocystis.Lys_pathway.txt \
examples/Synechocystis.model.tab \
examples/Synechocystis.model_drgs.tab \
results/Synechocystis_Lys_nem.csv
```

#### Example: MDF analysis on TCA cycle

```
./mdf.py --min_conc 0.0000001 --max_conc 0.1 \
--constraints examples/tca.concentrations.tab \
--ratios examples/tca.ratios.tab examples/tca.model.tab \
examples/tca.model_drgs.tab results/tca_mdf.csv
```

<a name="hit"></a>
## Hit-and-run analysis

`sampling.py` performs hit-and-run sampling of feasible metabolite concentration sets by a random walk through the solution space. Run `./sampling.py -h` to list all options of the script.

#### Quick example: Hit-and-run sampling on TCA cycle starting from and maintaining MDF

```
./sampling.py --reactions examples/tca.model.tab \
--std_drG examples/tca.model_drgs.tab \
--outfile results/tca.mdf_sampling.tab \
--constraints examples/tca.concentrations.tab \
--concs results/tca_mdf.csv --steps 1000 --mdf 1.3
```

#### Advanced example: Hit-and-run sampling on TCA cycle with replicates

Here we sample TCA cycle metabolite concentrations five times, computed in parallel and always starting from the MDF, saving every 10th step. In the first sampling all reactions must have a driving force > 0 kJ/mol, while in the second sampling all reactions must have a driving force > 1.18 kJ/mol, which is 90% of the MDF.

```
# Sample 5 replicates of 10000 steps in parallel at feasible delta G
for i in {1..5}; do
  echo "./sampling.py --reactions examples/tca.model.tab \
  --std_drG examples/tca.model_drgs.tab \
  --ratios examples/tca.ratio_range.tab \
  --outfile results/tca_sampling/tca.Feasible.${i}.tab \
  --constraints examples/tca.concentrations.tab \
  --concs results/tca_mdf.csv --steps 10000 --max_conc 1.2 -n 10"
done | parallel --no-notice --jobs 5

# Sample 5 replicates of 10000 steps in parallel at optimum delta G (90% of MDF)
for i in {1..5}; do
  echo "./sampling.py --reactions examples/tca.model.tab \
  --std_drG examples/tca.model_drgs.tab \
  --ratios examples/tca.ratio_range.tab \
  --outfile results/tca_sampling/tca.Optimal.${i}.tab \
  --constraints examples/tca.concentrations.tab \
  --concs results/tca_mdf.csv --steps 10000 --max_conc 1.2 \
  --mdf 1.18 -n 10"
done | parallel --no-notice --jobs 5
```

##### Convert reaction table to stoichiometric matrix

`stoich.py` converts a model in reaction table format to a model in stoichiometric matrix format. Run `./stoich.py -h` to list all options of the script. This representation of the model is needed for plotting in the next step.

```
./stoich.py examples/tca.model.tab results/tca.stoich.tab
```

##### Plot the output of multiple hit-and-run sampling replicates

`plot_samples.R` creates principal component analysis (PCA), concentration, and driving force plots from hit-and-run sampling experiments. Run `./plot_samples.R -h` to list all options of the script. All sampling results files used as input must be stored without any other files in one directory (`--indir` option). Sampling results files must have the `<path>/<indir>/<experiment>.<group>.<replicate>.tab` name format (even if there is only one group and/or only one replicate).

```
./plot_samples.R -i results/tca_sampling -S results/tca.stoich.tab \
-G examples/tca.model_drgs.tab -c examples/tca.concentrations.tab \
-x C00001 -o results/tca_sampling_plots
```

Plotting yields the following output files:
```
results/tca_sampling_plots.sampling_pca.png
results/tca_sampling_plots.sampling_pca_combo.png
results/tca_sampling_plots.sampling_concs.pdf
results/tca_sampling_plots.sampling_concs_combo.pdf
results/tca_sampling_plots.sampling_dfs.pdf
results/tca_sampling_plots.sampling_dfs_combo.pdf
results/tca_sampling_plots.concs.pca.pdf
results/tca_sampling_plots.dfs.pca.pdf
```

The `combo` tag indicates that replicates have been combined in the plot and are not displayed individually, as in the other variant of each plot type. Combining multiple  replicates improves coverage of the solution space.

| `results/tca_sampling_plots.sampling_pca.png` |
| --- |
| ![alt text](examples/images/tca_example_pca.png "Concentration PCA visualization of random walk for one replicate") |
| The PCA divided by individual samples shows how each random walk in the hit-and-run algorithm has traversed the solution space. The point color indicates the index of each feasible metabolite concentration set (fMCS) produced by the random walk. The figure below shows a cutout with a single replicate, where the left panel is the random walk at a driving force > 0 kJ/mol ("Feasible"), and the right panel has a driving force > 1.18 kJ/mol ("Optimal"). Note that the actual output file contains all replicates and is provided in the `results` directory. |

| `results/tca_sampling_plots.sampling_pca_combo.png` |
| --- |
| ![alt text](examples/images/tca_example_pca_combo.png "Concentration PCA visualization of random walk for all five replicates") |
| The combined PCA plot layers all replicate random walks, giving an indication of the coverage of the solution space. Each replicate up to 12 in total is given a unique color. The transparency attempts to show the step index. |

| `results/tca_sampling_plots.sampling_concs.pdf` |
| --- |
| ![alt text](examples/images/tca_example_concs.png "Visualization of random walk concentration distributions for all five replicates") |
| Concentration distributions obtained from the random walk are first visualized divided by replicate. "Feasible" is a driving force > 0 kJ/mol and "Optimal" is a driving force > 1.18 kJ/mol. Note that "Feasible" and "Optimal" are the `Group` variable names we gave to the two experiments through their file names. Up to 9 such names/groups can be selected freely to accomodate different types of experiments. In the concentration plots below, the vertical line marks on the x axes indicate the boundaries of the input concentration ranges. |

| `results/tca_sampling_plots.sampling_concs_combo.pdf` |
| --- |
| ![alt text](examples/images/tca_example_concs_combo.png "Visualization of random walk concentration distributions for all replicates combined") |
| As the random walks are quite noisy, it is best to combine multiple replicates to obtain smoother distributions. It can also be helpful to perform more steps, for example in the millions. The combined plot visualizes all values from all replicates of each group. Distributions are compared between `Group` variables with the Kolmogorov-Smirnov algorithm (function `ks.test`), providing the _D_ statistic that measures how different two distributions are. The _D_ statistic is indicated by the color ("Difference") of vertical lines connecting two distributions, situated in each facet to the right of the distributions. Any two distributions with _D_ ≤ 0.15 are considered equal and therefore do not get a connecting vertical line. |

| `results/tca_sampling_plots.sampling_dfs.pdf` |
| --- |
| ![alt text](examples/images/tca_example_dfs.png "Visualization of random walk driving force distributions for all five replicates") |
| Finally, the plotting script calculates the driving forces for each reaction based on the sampled metabolite concentrations. The driving forces are first plotted per replicate. Note that the x axis uses a square root scale. |

| `results/tca_sampling_plots.sampling_dfs_combo.pdf` |
| --- |
| ![alt text](examples/images/tca_example_dfs_combo.png "Visualization of random walk driving force distributions for all replicates combined") |
| Then, driving forces are plotted with replicates combined to yield a smoother representation of the thermodynamic solution space. It appears that [R01082](https://www.genome.jp/dbget-bin/www_bget?rn:R01082), _i.e._ fumarate hydratase, is somewhat of a thermodynamic bottleneck in this experiment. As for the concentrations, the Kolmogorov-Smirnov _D_ statistic ("Difference") is used to highlight differences in the distribution shapes. |

| `results/tca_sampling_plots.concs.pca.pdf` | `results/tca_sampling_plots.dfs.pca.pdf` |
| --- | --- |
| ![alt text](examples/images/tca_example_concs_pca.png "Visualization of metabolite co-variation through PCA of all concentrations") | ![alt text](examples/images/tca_example_dfs_pca.png "Visualization of reaction co-variation through PCA of all driving forces") |
| PCA of all metabolite concentrations reveals co-variation between metabolites. | PCA of all driving forces reveals co-variation between reactions. |

##### Use a plot config file to filter, order, and rename reactions, metabolites and groups

The `plot_samples.R` script can be supplied with a config file to filter, order, and rename reactions and metabolites. The example config file `examples/tca.plot_config.tab` is provided. Let's look at the contents:

```
cat examples/tca.plot_config.tab | column -tn -s $'\t'
```

```
Id        Name                      Type
R00351    Citrate synthase          Reaction
R01325    Citrate hydro-lyase       Reaction
R01900    Isocitrate hydro-lyase    Reaction
R00709    Isocitrate dehydrogenase  Reaction
R08549    2OG dehydrogenase         Reaction
R00405    Succinyl-CoA synthetase   Reaction
M00148    Succinate dehydrogenase   Reaction
R01082    Fumarate hydratase        Reaction
R00342    Malate dehydrogenase      Reaction
C00122    Fumarate                  Metabolite
C00149    Malate                    Metabolite
C00311    Isocitrate                Metabolite
C00036    Oxaloacetate              Metabolite
C00042    Succinate                 Metabolite
Feasible  Free                      Group
Optimal   MDF                       Group
```

The config file has three columns:
- `Id` is the identifier used for reactions, metabolites, or groups in the model and results files.
- `Name` is the alternative name we want to display in the plots, replacing the identifiers that are otherwise displayed.
- `Type` is one of `Reaction`, `Metabolite`, or `Group`, and defines what variable the identifier and name refers to. At least one line for each of the three types must be present in the config file to create meaningful plots and avoid crashes.

The exact order of reactions, metabolites, and groups in the config file is used to order the plots in the output accordingly. Thus, we use the plot config file to get nicely ordered TCA cycle plots with meaningful names:

```
./plot_samples.R -i results/tca_sampling -S results/tca.stoich.tab \
-G examples/tca.model_drgs.tab -c examples/tca.concentrations.tab \
-C examples/tca.plot_config.tab -o results/tca_cfg/tca_cfg_sampling_plots
```

These are the plots made with the help of the config file:

```
results/tca_cfg/tca_cfg_sampling_plots.sampling_pca.png
results/tca_cfg/tca_cfg_sampling_plots.sampling_pca_combo.png
results/tca_cfg/tca_cfg_sampling_plots.sampling_concs.pdf
results/tca_cfg/tca_cfg_sampling_plots.sampling_concs_combo.pdf
results/tca_cfg/tca_cfg_sampling_plots.sampling_dfs.pdf
results/tca_cfg/tca_cfg_sampling_plots.sampling_dfs_combo.pdf
results/tca_cfg/tca_cfg_sampling_plots.concs.pca.pdf
results/tca_cfg/tca_cfg_sampling_plots.dfs.pca.pdf
```

The concentrations now display only the five selected metabolites, C00122, C00149, C00311, C00036, and C00042, using their names Fumarate, Malate, Isocitrate, Oxaloacetate, and Succinate. Note that the groups of samples have been renamed to Free and MDF.

![alt text](examples/images/tca_cfg_example_concs.png "Re-configured concentration distribution visualization")

The driving forces are now displayed under more informative reaction names ordered according to their position in the TCA cycle.

![alt text](examples/images/tca_cfg_example_dfs.png "Re-configured driving force distributino visualization")

The metabolite concentrations and reaction driving forces PCA plots show the data of interest according to the definitions in the plot config file.

| ![alt text](examples/images/tca_cfg_example_concs_pca.png "Re-configured concentrations PCA plot") | ![alt text](examples/images/tca_cfg_example_dfs_pca.png "Re-configured driving forces PCA plot") |
| --- | --- |

<a name="syntrophism"></a>
## Case study: Syntrophic communities

One application of _thermosampler_ is the study of the metabolite concentrations and thermodynamics of metabolic reactions in microbial communities. Syntrophic propionate oxidation and methanogenesis was considered in this case study.

### Syntrophic propionate oxidizer model

A model was prepared for a syntrophic propionate oxidizing bacterium (SPOB) based on [_Candidatus_ Syntrophopropionicum ammoniitolerans](https://doi.org/10.1111/1462-2920.15388) (**[Fig 2](#fig2)**). The model was parameterized using data from literature and [Equilibrator](https://equilibrator.weizmann.ac.il/).

| Model component | File |
| --- | --- |
| Reactions describing metabolic network | `data/methylmalonyl.model.tab` |
| Stoichiometric matrix from `stoich.py` | `data/methylmalonyl.stoich.tab` |
| Standard reaction Gibbs free energy changes | `data/methylmalonyl.model_drGs.tab` |
| Default metabolite concentration ranges | `data/methylmalonyl.concentrations.tab` |
| Fixed concentration ratios for MDF | `data/methylmalonyl.ratios.tab` |
| Concentration ratio ranges for sampling | `data/methylmalonyl.ratio_range.tab` |

<a name="fig2"></a>

| ![alt text](examples/images/methylmalonyl.model.png "Model of SPOB") |
| --- |
| **Fig 2. Model of syntrophic propionate oxidizing bacterium.** Protons (H<sup>+</sup>) are not shown. Fd<sub>red</sub>, reduced ferredoxin; Fd<sub>ox</sub>, oxidized ferredoxin; MK, menaquinone. |

Fixed extracellular high and low propionate (input metabolite) and acetate (output metabolite) concentrations were defined to allow _thermosampler_ analysis regarding the thermodynamic effects of such conditions.

| Concentrations | File |
| --- | --- |
| Default | `data/mmcoa_fixed_prp_ac/methylmalonyl.concentrations.free.tab` |
| High acetate | `data/mmcoa_fixed_prp_ac/methylmalonyl.concentrations.acHi.tab` |
| Low propionate | `data/mmcoa_fixed_prp_ac/methylmalonyl.concentrations.prpLo.tab` |
| High prp/Low ac | `data/mmcoa_fixed_prp_ac/methylmalonyl.concentrations.prpHi_acLo.tab` |
| Low prp/High ac | `data/mmcoa_fixed_prp_ac/methylmalonyl.concentrations.prpLo_acHi.tab` |

The SPOB model was used to sample thermodynamically feasible metabolite concentrations for all the conditions above as described in this Bash script:

```
source/methylmalonyl_fixed_prp_ac.sh
```

The results were plotted using `plot_samples.R` and yielded the following output:

<details>
<summary>Full dataset plots.</summary>

```
# Random walk PCA
results/mmcoa_fixed_prp_ac_plots.sampling_pca_combo.png
results/mmcoa_fixed_prp_ac_plots.sampling_pca.png

# Concentrations and driving forces PCA
results/mmcoa_fixed_prp_ac_plots.concs.pca.pdf
results/mmcoa_fixed_prp_ac_plots.dfs.pca.pdf

# Concentration distributions
results/mmcoa_fixed_prp_ac_plots.sampling_concs_combo.pdf
results/mmcoa_fixed_prp_ac_plots.sampling_concs.pdf

# Driving force distributions
results/mmcoa_fixed_prp_ac_plots.sampling_dfs_combo.pdf
results/mmcoa_fixed_prp_ac_plots.sampling_dfs.pdf
```
</details>

<details>
<summary>Reduced dataset plots.</summary>

```
# Random walk PCA
results/mmcoa_fixed_prp_ac_plots_reduced.sampling_pca_combo.png
results/mmcoa_fixed_prp_ac_plots_reduced.sampling_pca.png

# Concentrations and driving forces PCA
results/mmcoa_fixed_prp_ac_plots_reduced.concs.pca.pdf
results/mmcoa_fixed_prp_ac_plots_reduced.dfs.pca.pdf

# Concentration distributions
results/mmcoa_fixed_prp_ac_plots_reduced.sampling_concs_combo.pdf
results/mmcoa_fixed_prp_ac_plots_reduced.sampling_concs.pdf

# Driving force distributions
results/mmcoa_fixed_prp_ac_plots_reduced.sampling_dfs_combo.pdf
results/mmcoa_fixed_prp_ac_plots_reduced.sampling_dfs.pdf
```
</details>

### Syntrophic methanogenesis model

The SPOB was given company by an acetoclastic methanogen and a hydrogenotrophic methanogen to produce a multi-organism syntrophic model. To distinguish the different organisms' reactions, a prefix was introduced; `po_` for the SPOB, and `ac_` and `hm_` for the acetoclastic and hydrogenotrophic methanogens. The three organisms affected eachother through the overlap in metabolites, _i.e._ hydrogen, formate, acetate, methane, and carbon dioxide. The model was parameterized using data from literature and [Equilibrator](https://equilibrator.weizmann.ac.il/).

| Model component | File |
| --- | --- |
| Reactions describing metabolic network | `data/syntrophic.model.tab` |
| Stoichiometric matrix from `stoich.py` | `data/syntrophic.stoich.tab` |
| Standard reaction Gibbs free energy changes | `data/syntrophic.model_drGs.tab` |
| Default metabolite concentration ranges | `data/syntrophic.concentrations.tab` |
| Individual concentration sums for each organism | `data/syntrophic.conc_sums.tab` |
| Fixed concentration ratios for MDF | `data/syntrophic.ratios.tab` |
| Concentration ratio ranges for sampling | `data/syntrophic.ratio_range.tab` |

Additionally, a file with a fixed high methane concentration was prepared to make contrast with free methane concentration thermodynamics:

```
data/syntrophic.concentrations.high_CH4.tab
```

The syntrophic model was used to sample thermodynamically feasible metabolite concentrations in free and high methane concentration conditions, as described in these Bash scripts:

```
source/syntrophic_2x10x10M.sh
source/syntrophic_2x10x10M.high_CH4.sh
```

The results were plotted using `plot_samples.R` and yielded the following output:

<details open>
<summary>Methane comparison plots.</summary>

```
# Random walk PCA
results/syntrophic_CH4.sampling_pca_combo.png
results/syntrophic_CH4.sampling_pca.png

# Concentrations and driving forces PCA
results/syntrophic_CH4.concs.pca.pdf
results/syntrophic_CH4.dfs.pca.pdf

# Concentration distributions
results/syntrophic_CH4.sampling_concs_combo.pdf
results/syntrophic_CH4.sampling_concs.pdf

# Driving force distributions
results/syntrophic_CH4.sampling_dfs_combo.pdf
results/syntrophic_CH4.sampling_dfs.pdf
```
</details>


<a name="author"></a>
## Author

Johannes Asplund-Samuelsson (johannes.aspsam@gmail.com)
