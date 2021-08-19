#!/usr/bin/env Rscript

### COMMAND LINE ARGUMENTS #####################################################

library(optparse)

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
outprefix = opt$outprefix

# Required files must be supplied
if (sum(is.na(c(indir, stoich_file, drgs_file, concs_file, outprefix)))) {
  stop("Missing required arguments. See usage information (--help).")
}

# Load libraries
library(tidyverse)
library(foreach)
library(doMC)

# Load data
registerDoMC(detectCores())

exclude = str_split(exclude_string, ",") %>% unlist()

sampling = bind_rows(
  foreach(
    f=list.files(indir, full.names=T)
  ) %dopar% {
    Grp = str_split(f, "\\.")[[1]][2]
    Rep = str_split(f, "\\.")[[1]][3]
    read_tsv(f) %>% mutate(Group = Grp, Replicate = Rep)
  }
)

S = read_tsv(stoich_file) %>%
  rename(Metabolite = X1)
G = read_tsv(drgs_file, col_names=c("Reaction", "drG"))
RT = 8.31e-3 * 298.15

sample_groups = unique(sampling$Group)

# Numbers for formatting
n_rep = length(unique(sampling$Replicate)) # Number of replicates
n_grp = length(sample_groups) # Number of groups
n_met = nrow(S) - sum(!is.na(exclude)) # Number of metabolites
n_rxn = ncol(S) - 1

conc_ranges = concs_file %>%
  read_tsv(col_names=c("Metabolite", "Low", "High")) %>%
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
  filter(!(Metabolite %in% exclude))

# Perform PCA to plot random walk steps through solution space
sampling_X = sampling %>%
  select(-Group, -Run, -fMCS, -Replicate) %>%
  as.matrix()

sampling_pca = prcomp(scale(sampling_X, scale=F))

# Calculate fraction of variance per PC
library(scales)
sampling_pca_var = percent(
  sampling_pca$sdev^2 / sum(sampling_pca$sdev^2),
  accuracy = 0.1
)

# Create plotting dataframe
sampling_pca_plot = as_tibble(sampling_pca$x) %>%
  select(PC1, PC2) %>%
  bind_cols(select(sampling, Group, fMCS, Replicate))

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
gp = gp + scale_colour_brewer(palette="Paired", guide=F)
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
  mutate(Concentration = exp(lnC)*1000)

library(ggridges)

gp = ggplot(
  filter(concs, !(Metabolite %in% exclude)),
  aes(x=Concentration, fill=Group, y=paste(Group, Replicate))
)
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Metabolite~., ncol=5)
gp = gp + scale_fill_brewer(palette="YlGnBu", guide=F)
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

gp = ggplot(
  filter(concs, !(Metabolite %in% exclude)),
  aes(x=Concentration, fill=Group, y=Group)
)
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Metabolite~., ncol=5)
gp = gp + scale_fill_brewer(palette="YlGnBu", guide=F)
gp = gp + geom_point(
  data=conc_ranges, shape=124, fill="black", color="black"
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
  paste(outprefix, ".sampling_concs_combo.pdf", sep=""),
  width=210/25.4,
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

dfs = sampling_G %>% gather(Reaction, DF, -Group, -fMCS, -Replicate)

gp = ggplot(dfs, aes(x=DF, fill=Group, y=paste(Group, Replicate)))
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Reaction~., ncol=5)
gp = gp + scale_fill_brewer(palette="YlGnBu", guide=F)
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

gp = ggplot(dfs, aes(x=DF, fill=Group, y=Group))
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Reaction~., ncol=5)
gp = gp + scale_fill_brewer(palette="YlGnBu", guide=F)
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
  width=210/25.4,
  height=230/25.4/(44/5*2)*n_rxn/5*n_grp + 50/25.4,
  limitsize = FALSE
)
