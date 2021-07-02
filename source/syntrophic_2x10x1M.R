options(width=110)
library(tidyverse)
library(foreach)
library(doMC)

# Load data
registerDoMC(20)
sampling = bind_rows(
  foreach(
    f=list.files("results/syntrophic_2x10x1M", full.names=T)
  ) %dopar% {
    # Read only every 10 rows
    fp = pipe(paste("(head -1", f, "; sed -n '2~10p'", f, ")"))
    Lim = ifelse(grepl("mdf", f), "Optimum", "Feasible")
    Rep = str_split(f, "\\.")[[1]][3]
    read_tsv(fp) %>% mutate(Limit = Lim, Replicate = Rep)
  }
)

S = read_tsv("data/syntrophic.stoich.tab") %>%
  rename(Metabolite = X1) %>%
  filter(Metabolite != "H2O")
G = read_tsv(
  "data/syntrophic.model_drGs.tab",
  col_names=c("Reaction", "drG")
)
RT = 8.31e-3 * 298.15

conc_ranges = read_tsv(
  "data/syntrophic.concentrations.tab",
  col_names=c("Metabolite", "Low", "High")
) %>%
  gather(Bound, Concentration, -Metabolite) %>%
  mutate(Concentration = 1000*Concentration) %>%
  inner_join(
    tibble(
      Limit=rep(c("Feasible", "Optimum"), each=2),
      Bound=rep(c("Low", "High"), 2)
    )
  ) %>%
  filter(Metabolite != "H2O")

# Perform PCA to plot random walk steps through solution space
sampling_X = sampling %>%
  select(-Limit, -Run, -fMCS, -Replicate, -H2O) %>%
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

ggsave("results/syntrophic.sampling_pca.2x10x1M.png", w=7, h=27, dpi=250)

sampling_pca_plot_small = filter(sampling_pca_plot, fMCS %% 100 == 0)

gp = ggplot(
  sampling_pca_plot,
  aes(x=PC1, y=PC2, alpha=fMCS, colour=Replicate)
)
gp = gp + geom_point(size=0.1, stroke=0)
gp = gp + scale_alpha_continuous(range=c(0.5,0.1))
gp = gp + scale_colour_manual(
  values=c(
    "#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3",
    "#003c30","#01665e","#35978f","#80cdc1","#c7eae5"),
  guide=F
)
gp = gp + labs(
            x=paste("PC1 (", sampling_pca_var[1], ")", sep=""),
            y=paste("PC2 (", sampling_pca_var[2], ")", sep="")
          )
gp = gp + theme_bw()
gp = gp + facet_grid(.~Limit)
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  strip.background = element_blank(),
  aspect.ratio = 1
)

ggsave("results/syntrophic.sampling_pca_compact.2x10x1M.png", w=7, h=4, dpi=250)

# Plot concentration ranges
concs = sampling %>%
  select(-Run) %>%
  gather(Metabolite, lnC, -Limit, -fMCS, -Replicate) %>%
  mutate(Concentration = exp(lnC)*1000) %>%
  filter(Metabolite != "H2O")

library(ggridges)

small_concs = filter(concs, fMCS %% 100 == 0)

gp = ggplot(
  #small_concs,
  concs,
  aes(
    x=Concentration,
    fill=Limit,
    y=paste(Limit, Replicate)
  )
)
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Metabolite~., ncol=5)
gp = gp + scale_fill_manual(values=c("#e08214","#8073ac"), guide=F)
gp = gp + geom_point(
  data = conc_ranges %>% inner_join(
    tibble(
      Replicate=rep(1:10, each=2),
      Bound=rep(c("Low", "High"), 10)
    )
  ),
  shape=124, fill="black", color="black"
)
gp = gp + scale_x_log10()
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
  "results/syntrophic.sampling_concs.2x10x1M.pdf", gp, w=210/25.4, h=1200/25.4
)

gp = ggplot(
  #small_concs,
  concs,
  aes(
    x=Concentration,
    fill=Limit,
    y=Limit
  )
)
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Metabolite~., ncol=5)
gp = gp + scale_fill_manual(values=c("#e08214","#8073ac"), guide=F)
gp = gp + geom_point(
  data=conc_ranges, shape=124, fill="black", color="black"
)
gp = gp + scale_x_log10()
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
  "results/syntrophic.sampling_concs_combo.2x10x1M.pdf",
  gp, w=210/25.4, h=420/25.4
)

# Calculate driving forces

# Arrange the stoichiometric matrix correctly
S = S %>% arrange(match(Metabolite, colnames(sampling_X)))
G = G %>% arrange(match(Reaction, colnames(select(S, -Metabolite))))

sampling_G = as_tibble(
  -1*(sampling_X %*% as.matrix(select(S, -Metabolite)) * RT +
  outer(rep.int(1L, nrow(sampling_X)), G$drG))
) %>% bind_cols(select(sampling, Limit, fMCS, Replicate))

dfs = sampling_G %>% gather(Reaction, DF, -Limit, -fMCS, -Replicate)

small_dfs = filter(dfs, fMCS %% 100 == 0)

gp = ggplot(
  #small_dfs,
  dfs,
  aes(
    x=DF,
    fill=Limit,
    y=paste(Limit, Replicate)
  )
)
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Reaction~., ncol=5)
gp = gp + scale_fill_manual(values=c("#e08214","#8073ac"), guide=F)
gp = gp + scale_x_sqrt()
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
  "results/syntrophic.sampling_dfs.2x10x1M.pdf", gp, w=210/25.4, h=800/25.4
)

gp = ggplot(
  #small_dfs,
  dfs,
  aes(
    x=DF,
    fill=Limit,
    y=Limit
  )
)
gp = gp + geom_density_ridges(alpha=0.5)
gp = gp + facet_wrap(Reaction~., ncol=5)
gp = gp + scale_fill_manual(values=c("#e08214","#8073ac"), guide=F)
gp = gp + scale_x_sqrt()
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
  "results/syntrophic.sampling_dfs_combo.2x10x1M.pdf",
  gp, w=210/25.4, h=280/25.4
)
