options(width=110)
library(tidyverse)

# Load data
sampling = bind_rows(
  # Replicate 1
  read_tsv("results/methylmalonyl.sampling.tab") %>%
    mutate(Limit = "Feasible", Replicate = 1),
  read_tsv("results/methylmalonyl.mdf_sampling.tab") %>%
    mutate(Limit = "Optimum", Replicate = 1),
  # Replicate 2
  read_tsv("results/methylmalonyl.sampling_rep2.tab") %>%
    mutate(Limit = "Feasible", Replicate = 2),
  read_tsv("results/methylmalonyl.mdf_sampling_rep2.tab") %>%
    mutate(Limit = "Optimum", Replicate = 2),
  # Replicate 3
  read_tsv("results/methylmalonyl.sampling_rep3.tab") %>%
    mutate(Limit = "Feasible", Replicate = 3),
  read_tsv("results/methylmalonyl.mdf_sampling_rep3.tab") %>%
    mutate(Limit = "Optimum", Replicate = 3),
) %>% select(Limit, everything())

S = read_tsv("data/methylmalonyl.stoich.tab") %>% rename(Reaction = X1)

# Perform PCA to plot random walk steps through solution space
sampling_X = sampling %>%
  select(-Limit, -Run, -fMCS, -Replicate) %>%
  as.matrix()

# Remove constant columns
sampling_X = sampling_X[,apply(sampling_X, 2, var) != 0]

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
  bind_cols(select(sampling, Limit, fMCS, Replicate)) %>%
  mutate(Limit = factor(Limit, levels=c("Optimum", "Feasible")))

gp = ggplot(
  sampling_pca_plot,
  aes(x=PC1, y=PC2, colour=fMCS)
)
gp = gp + geom_point(alpha=0.2, size=0.5, stroke=0)
gp = gp + scale_colour_distiller(palette="PuBuGn")
gp = gp + labs(
            x=paste("PC1 (", sampling_pca_var[1], ")", sep=""),
            y=paste("PC2 (", sampling_pca_var[2], ")", sep="")
          )
gp = gp + theme_bw()
gp = gp + facet_grid(Replicate~Limit)
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  aspect.ratio = 1
)

ggsave("results/methylmalonyl.sampling_pca.png", w=7, h=9.5, dpi=250)

# Plot concentration ranges
concs = sampling %>%
  select(-Run) %>%
  gather(Metabolite, lnC, -Limit, -fMCS, -Replicate) %>%
  mutate(Concentration = exp(lnC)*1000) %>%
  filter(Metabolite != "H2O")

library(ggridges)

gp = ggplot(
  #filter(concs, fMCS %% 10 == 0),
  concs,
  aes(
    x=Concentration,
    fill=paste(Limit, Replicate),
    y=paste(Limit, Replicate)
  )
)
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Metabolite~., ncol=5)
gp = gp + scale_fill_brewer(palette="PuOr", guide=F)
gp = gp + scale_x_log10()
gp = gp + theme_bw()
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  axis.text.x = element_text(angle=60, vjust=1, hjust=1),
  axis.title.y = element_blank()
)

ggsave("results/methylmalonyl.sampling_concs.pdf", gp, w=210/25.4, h=297/25.4)
