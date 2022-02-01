#!/usr/bin/env Rscript

# Load libraries
library(tidyverse)
library(foreach)
library(doMC)
library(ggridges)
library(scales)
library(optparse)
library(ggrepel)

### COMMAND LINE ARGUMENTS #####################################################

# List command line options
option_list = list(
  make_option(
    c("-i", "--indir"), type="character", default=NA,
    help="Sampling results directory (required)."
  ),
  make_option(
    c("-S", "--stoich"), type="character", default=NA,
    help="Stoichiometric matrix file (required)."
  ),
  make_option(
    c("-G", "--drgs"), type="character", default=NA,
    help="Standard Gibbs free energy change value (delta G) file (required)."
  ),
  make_option(
    c("-c", "--concs"), type="character", default=NA,
    help="Concentration range file (required)."
  ),
  make_option(
    c("-x", "--exclude"), type="character", default="",
    help="Metabolites to exclude from plotting, separated by comma."
  ),
  make_option(
    c("-r", "--remove"), type="character", default="",
    help="Reactions to remove from plotting, separated by comma."
  ),
  make_option(
    c("-C", "--config"), type="character", default="",
    help="Table to rename, filter, and order Group, Reaction, and Metabolite."
  ),
  make_option(
    c("-o", "--outprefix"), type="character", default=NA,
    help="Prefix for output (required)."
  )
)

# Parse command line arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Save options into new variables
indir = opt$indir
stoich_file = opt$stoich
drgs_file = opt$drgs
concs_file = opt$concs
exclude_string = opt$exclude
remove_string = opt$remove
config_file = opt$config
outprefix = opt$outprefix

# Required files must be supplied
if (sum(is.na(c(indir, stoich_file, drgs_file, concs_file, outprefix)))) {
  stop("Missing required arguments. See usage information (--help).")
}

# Load data
registerDoMC(detectCores())


sampling = bind_rows(
  foreach(
    f=list.files(indir, full.names=T)
  ) %dopar% {
    Grp = str_split(f, "\\.")[[1]][2]
    Rep = str_split(f, "\\.")[[1]][3]
    read_tsv(f) %>% mutate(Group = Grp, Replicate = Rep)
  }
)

S = read_tsv(stoich_file) %>% rename(Metabolite = 1)
G = read_tsv(drgs_file, comment='#', col_names=c("Reaction", "drG"))
RT = 8.31e-3 * 298.15


# Load configuration data frame and exclude metabolites and reactions
if (config_file != "") {
  config_df = read_tsv(config_file, comment='#')
} else {
  config_df = tibble(
    Id = c(G$Reaction, S$Metabolite, unique(sampling$Group)),
    Name = c(G$Reaction, S$Metabolite, unique(sampling$Group)),
    Type = c(
      rep("Reaction", length(G$Reaction)),
      rep("Metabolite", length(S$Metabolite)),
      rep("Group", length(unique(sampling$Group)))
    )
  )
}

exclude_metabolites = str_split(exclude_string, ",") %>% unlist()
exclude_reactions = str_split(remove_string, ",") %>% unlist()

config_df = config_df %>%
  filter(!(Type == "Metabolite" & Id %in% exclude_metabolites)) %>%
  filter(!(Type == "Reaction" & Id %in% exclude_reactions)) %>%
  mutate(Name = factor(Name, levels=Name))

sample_groups = config_df %>% filter(Type == "Group") %>% pull(Id)

# Numbers for formatting
n_rep = length(unique(sampling$Replicate)) # Number of replicates
n_grp = length(sample_groups) # Number of groups
n_met = nrow(filter(config_df, Type == "Metabolite")) # Number of metabolites
n_rxn = nrow(filter(config_df, Type == "Reaction")) # Number of reactions

conc_ranges = concs_file %>%
  read_tsv(col_names=c("Metabolite", "Low", "High"), comment='#') %>%
  gather(Bound, Concentration, -Metabolite) %>%
  mutate(Concentration = 1000*Concentration) %>%
  # Add every combination of Group and Bound
  inner_join(
    tibble(
      Group=rep(sample_groups, each=2),
      Bound=rep(c("Low", "High"), n_grp)
    )
  ) %>%
  # Remove metabolites that are excluded
  inner_join(
    config_df %>%
      filter(Type == "Metabolite") %>%
      select(Id, Name) %>%
      rename(Metabolite = Id)
  ) %>%
  select(-Metabolite) %>%
  rename(Metabolite = Name) %>%
  inner_join(
    config_df %>%
      filter(Type == "Group") %>%
      select(Id, Name) %>%
      rename(Group = Id)
  ) %>%
  select(-Group) %>%
  rename(Group = Name)

# Perform PCA to plot random walk steps through solution space
sampling_X = sampling %>%
  select(-Group, -Run, -fMCS, -Replicate) %>%
  as.matrix()

sampling_pca = prcomp(
  scale(sampling_X[,apply(sampling_X, 2, var, na.rm=TRUE) != 0], scale=T)
)

# Calculate fraction of variance per PC
sampling_pca_var = percent(
  sampling_pca$sdev^2 / sum(sampling_pca$sdev^2),
  accuracy = 0.1
)

# Create plotting dataframe
sampling_pca_plot = as_tibble(sampling_pca$x) %>%
  select(PC1, PC2) %>%
  bind_cols(select(sampling, Group, fMCS, Replicate))

# Rename and subset Group according to config
sampling_pca_plot = filter(config_df, Type == "Group") %>%
  select(-Type) %>%
  rename(Group=Id) %>%
  inner_join(sampling_pca_plot) %>%
  select(-Group) %>%
  rename(Group=Name)

gp = ggplot(
  sampling_pca_plot,
  aes(x=PC1, y=PC2, colour=fMCS)
)
gp = gp + geom_point(alpha=0.2, size=1, stroke=0)
gp = gp + scale_colour_distiller(palette="PuBuGn")
gp = gp + labs(
            x=paste("PC1 (", sampling_pca_var[1], ")", sep=""),
            y=paste("PC2 (", sampling_pca_var[2], ")", sep="")
          )
gp = gp + theme_bw()
gp = gp + facet_grid(Replicate~Group)
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  aspect.ratio = 1
)

ggsave(
  paste(outprefix, ".sampling_pca.png", sep=""),
  width=7/2*n_grp,
  height=27/10*n_rep,
  dpi=250,
  limitsize = FALSE
)

gp = ggplot(
  sampling_pca_plot,
  aes(x=PC1, y=PC2, alpha=fMCS, colour=Replicate)
)
gp = gp + geom_point(size=0.5, stroke=0)
gp = gp + scale_alpha_continuous(range=c(0.5,0.1))
gp = gp + scale_colour_brewer(palette="Paired", guide="none")
gp = gp + labs(
            x=paste("PC1 (", sampling_pca_var[1], ")", sep=""),
            y=paste("PC2 (", sampling_pca_var[2], ")", sep="")
          )
gp = gp + theme_bw()
gp = gp + facet_grid(.~Group)
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  aspect.ratio = 1
)

ggsave(
  paste(outprefix, ".sampling_pca_combo.png", sep=""),
  width=7/2*n_grp, height=4, dpi=250,
  limitsize = FALSE
)

# Plot concentration ranges
concs = sampling %>%
  select(-Run) %>%
  gather(Metabolite, lnC, -Group, -fMCS, -Replicate) %>%
  mutate(Concentration = exp(lnC)*1000) %>%
  inner_join(
    config_df %>%
      filter(Type == "Metabolite") %>%
      select(Id, Name) %>%
      rename(Metabolite = Id)
  ) %>%
  select(-Metabolite) %>%
  rename(Metabolite = Name) %>%
  inner_join(
    config_df %>%
      filter(Type == "Group") %>%
      select(Id, Name) %>%
      rename(Group = Id)
  ) %>%
  select(-Group) %>%
  rename(Group = Name)

gp = ggplot(
  concs,
  aes(
    x=Concentration, fill=Group,
    y=factor(
      paste(Group, Replicate),
      levels = paste(
        config_df %>%
          filter(Type == "Group") %>%
          pull(Name) %>%
          as.character() %>%
          rep(each=n_rep),
        concs %>%
          pull(Replicate) %>%
          unique()
      )
    )
  )
)
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Metabolite~., ncol=5)
gp = gp + scale_fill_brewer(palette="YlGnBu", guide="none")
gp = gp + geom_point(
  data = conc_ranges %>%
    inner_join(
      tibble(
        Replicate=rep(unique(sampling$Replicate), each=2),
        Bound=rep(c("Low", "High"), n_rep)
      )
    ),
  shape=124, fill="black", color="black"
)
gp = gp + scale_x_log10(breaks=10**(-5:3), minor_breaks=NULL)
gp = gp + theme_bw()
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  axis.text.x = element_text(angle=60, vjust=1, hjust=1),
  axis.title.y = element_blank()
)
gp = gp + xlab("Concentration (mM)")

ggsave(
  paste(outprefix, ".sampling_concs.pdf", sep=""),
  width=210/25.4,
  height=1150/25.4/(72/5*2*10)*n_met/5*n_rep*n_grp + 50/25.4,
  limitsize = FALSE
)

# Perform Kolmogorov-Smirnov test of distributions
test_pairs = concs$Group %>%
  unique() %>%
  as.character() %>%
  combn(2) %>%
  t() %>%
  as_tibble() %>%
  rename(A = V1, B = V2) %>%
  mutate(
    A = factor(A, levels=levels(concs$Group)),
    B = factor(B, levels=levels(concs$Group))
  )

compare_distributions = function(df, test_pairs){
  test_pairs$D = 0
  for (i in 1:nrow(test_pairs)){
    A = filter(df, Group == test_pairs[i,]$A)$Value
    B = filter(df, Group == test_pairs[i,]$B)$Value
    test_output = ks.test(A, B, alternative="two.sided")
    test_pairs[i,"D"] = test_output$statistic
  }
  return(test_pairs)
}

registerDoMC(4)

concs_comparison = bind_rows(
  foreach (x = unique(concs$Metabolite)) %dopar% {
    concs %>%
      filter(Metabolite == x) %>%
      select(Group, lnC) %>%
      rename(Value = lnC) %>%
      compare_distributions(test_pairs) %>%
      mutate(Metabolite = x)
  }
)

# Calculate where to place lines showing dissimilarity of distributions
min_conc = concs %>% pull(Concentration) %>% min()

max_conc = concs %>% pull(Concentration) %>% max()

max_min_conc = concs %>%
  group_by(Metabolite) %>%
  summarise(c_min=min(Concentration), c_max=max(Concentration)) %>%
  mutate(
    Space_below = log10(c_min) - log10(min_conc),
    Space_above = log10(max_conc) - log10(c_max)
  )

concs_comparison_plot = concs_comparison %>%
  inner_join(
    config_df %>%
      filter(Type == "Group") %>%
      select(Name) %>%
      rename(A=Name) %>%
      mutate(Y = rank(as.integer(A)))
  ) %>%
  inner_join(
    config_df %>%
      filter(Type == "Group") %>%
      select(Name) %>%
      rename(B=Name) %>%
      mutate(Y_end = rank(as.integer(B)))
  ) %>%
  left_join(max_min_conc) %>%
  group_by(Metabolite) %>%
  mutate(X_order = rank(Y + Y_end, ties.method="first")) %>%
  ungroup() %>%
  mutate(X = c_max * (10 ** (0.3*X_order + 0.5)), Group = A) %>%
  filter(D > 0.15) %>%
  rename(Difference = D)

# Plot combined distributions of concentrations
gp = ggplot(concs, aes(x=Concentration, fill=Group, y=Group))
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Metabolite~., ncol=5)
gp = gp + scale_fill_brewer(palette="YlGnBu", guide="none")
gp = gp + geom_point(
  data=conc_ranges, shape=124, fill="black", color="black"
)
gp = gp + geom_segment(
  data=concs_comparison_plot,
  mapping=aes(x=X, xend=X, y=Y, yend=Y_end, color = Difference),
)
gp = gp + scale_color_viridis_c(option="plasma")
gp = gp + scale_x_log10(
  breaks=10**(
    floor(log10(min_conc)):ceiling(log10(max(concs_comparison_plot$X)))
  ),
  minor_breaks=NULL
)
gp = gp + theme_bw()
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  axis.text.x = element_text(angle=60, vjust=1, hjust=1),
  axis.title.y = element_blank()
)
gp = gp + xlab("Concentration (mM)")

ggsave(
  paste(outprefix, ".sampling_concs_combo.pdf", sep=""),
  width=(210 + 15*nrow(test_pairs))/25.4,
  height=370/25.4/(72/5*2)*n_met/5*n_grp + 50/25.4,
  limitsize = FALSE
)

# Calculate driving forces

# Arrange the stoichiometric matrix correctly
S = S %>% arrange(match(Metabolite, colnames(sampling_X)))
G = G %>% arrange(match(Reaction, colnames(select(S, -Metabolite))))

sampling_G = as_tibble(
  -1*(sampling_X %*% as.matrix(select(S, -Metabolite)) * RT +
  outer(rep.int(1L, nrow(sampling_X)), G$drG))
) %>% bind_cols(select(sampling, Group, fMCS, Replicate))

# Organize driving forces with correct Group and Reaction selection
dfs = sampling_G %>%
  gather(Reaction, DF, -Group, -fMCS, -Replicate) %>%
  inner_join(
    config_df %>%
      filter(Type == "Reaction") %>%
      select(Id, Name) %>%
      rename(Reaction = Id)
  ) %>%
  select(-Reaction) %>%
  rename(Reaction = Name) %>%
  inner_join(
    config_df %>%
      filter(Type == "Group") %>%
      select(Id, Name) %>%
      rename(Group = Id)
  ) %>%
  select(-Group) %>%
  rename(Group = Name)

gp = ggplot(
  dfs,
  aes(
    x=DF, fill=Group,
    y=factor(
      paste(Group, Replicate),
      levels = paste(
        config_df %>%
          filter(Type == "Group") %>%
          pull(Name) %>%
          as.character() %>%
          rep(each=n_rep),
        concs %>%
          pull(Replicate) %>%
          unique()
      )
    )
  )
)
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Reaction~., ncol=5)
gp = gp + scale_fill_brewer(palette="YlGnBu", guide="none")
gp = gp + scale_x_sqrt(breaks=c(1,8,25,50,75,100), minor_breaks=NULL)
gp = gp + theme_bw()
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  axis.text.x = element_text(angle=60, vjust=1, hjust=1),
  axis.title.y = element_blank()
)
gp = gp + xlab("Driving force (kJ/mol)")

ggsave(
  paste(outprefix, ".sampling_dfs.pdf", sep=""),
  width=210/25.4,
  height=750/25.4/(44/5*2*10)*n_rxn/5*n_rep*n_grp + 50/25.4,
  limitsize = FALSE
)

# Compare distributions of DFs
dfs_comparison = bind_rows(
  foreach (x = unique(dfs$Reaction)) %dopar% {
    dfs %>%
      filter(Reaction == x) %>%
      select(Group, DF) %>%
      rename(Value = DF) %>%
      compare_distributions(test_pairs) %>%
      mutate(Reaction = x)
  }
)

# Calculate where to place lines showing dissimilarity of distributions
min_df = dfs %>% pull(DF) %>% min()

max_df = dfs %>% pull(DF) %>% max()

max_min_df = dfs %>%
  group_by(Reaction) %>%
  summarise(df_min=min(DF), df_max=max(DF))

dfs_comparison_plot = dfs_comparison %>%
  inner_join(
    config_df %>%
      filter(Type == "Group") %>%
      select(Name) %>%
      rename(A=Name) %>%
      mutate(Y = rank(as.integer(A)))
  ) %>%
  inner_join(
    config_df %>%
      filter(Type == "Group") %>%
      select(Name) %>%
      rename(B=Name) %>%
      mutate(Y_end = rank(as.integer(B)))
  ) %>%
  left_join(max_min_df) %>%
  group_by(Reaction) %>%
  mutate(X_order = rank(Y + Y_end, ties.method="first")) %>%
  ungroup() %>%
  mutate(
    X = (sqrt(df_max + max_df/20) + X_order * sqrt(max_df)/20) ** 2,
    Group = A
  ) %>%
  filter(!(Reaction %in% exclude_reactions)) %>%
  filter(D > 0.15) %>%
  rename(Difference = D)

gp = ggplot(dfs, aes(x=DF, fill=Group, y=Group))
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Reaction~., ncol=5)
gp = gp + scale_fill_brewer(palette="YlGnBu", guide="none")
gp = gp + geom_segment(
  data=dfs_comparison_plot,
  mapping=aes(x=X, xend=X, y=Y, yend=Y_end, color = Difference),
)
gp = gp + scale_color_viridis_c(option="plasma")
gp = gp + scale_x_sqrt(breaks=c(1,8,25,50,75,100), minor_breaks=NULL)
gp = gp + theme_bw()
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  axis.text.x = element_text(angle=60, vjust=1, hjust=1),
  axis.title.y = element_blank()
)
gp = gp + xlab("Driving force (kJ/mol)")

ggsave(
  paste(outprefix, ".sampling_dfs_combo.pdf", sep=""),
  width=(210 + 15*nrow(test_pairs))/25.4,
  height=230/25.4/(44/5*2)*n_rxn/5*n_grp + 50/25.4,
  limitsize = FALSE
)

# Plot PCA of metabolite concentrations and driving forces
plot_pca = function(pca_results, data_type){

  # Create plotting dataframes
  pca_plot = as.data.frame(pca_results$rotation)
  pca_plot$Id = rownames(pca_plot)

  # Select wanted data
  pca_plot = pca_plot %>%
    inner_join(filter(config_df, Type == data_type) %>% select(-Type))

  # Calculate fraction of variance per PC
  var_pc = percent(
    pca_results$sdev^2 / sum(pca_results$sdev^2),
    accuracy=0.1
  )

  gp = ggplot(pca_plot, aes(x=PC1, y=PC2, colour=PC3, label=Name))
  gp = gp + geom_point(aes(size=PC3))
  gp = gp + scale_size(range=c(1,3))
  gp = gp + scale_colour_gradient2(low="#d0d1e6", mid="#3690c0", high="#014636")
  gp = gp + geom_text_repel(force=3, size=4)
  gp = gp + labs(
              x=paste("PC1 (", var_pc[1], ")", sep=""),
              y=paste("PC2 (", var_pc[2], ")", sep=""),
              colour=paste("PC3 (", var_pc[3], ")", sep="")
            )
  gp = gp + theme_bw()

  tag = ifelse(data_type == "Reaction", "dfs", "concs")
  outfile = paste(outprefix, tag, "pca", "pdf", sep=".")

  ggsave(outfile, gp, height=15/2.54, width=18/2.54)

}

# Plot metabolites from concentration PCA
plot_pca(sampling_pca, "Metabolite")

# Create driving forces matrix
dfs_X = sampling_G %>%
  gather(Reaction, DF, -Group, -fMCS, -Replicate) %>%
  spread(Reaction, DF) %>%
  select(-fMCS, -Replicate, -Group) %>%
  as.matrix()

# Perform PCA on driving forces
df_pca = prcomp(scale(dfs_X, scale=T))

# Plot reactions from driving force PCA
plot_pca(df_pca, "Reaction")
