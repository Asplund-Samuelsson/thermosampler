options(width=110)
library(tidyverse)

# Load data
concs = read_tsv("results/toy.sampling.tab") %>%
  mutate(A = exp(A), B = exp(B))

start_end = concs[c(1, nrow(concs)),] %>%
  mutate(Label = c("Start", "Finish"))

gp = ggplot(concs, aes(x=A, y=B, color=fMCS))
gp = gp + geom_path(alpha=0.6)
gp = gp + geom_point(alpha=0.6)
gp = gp + geom_abline(slope=-1, intercept=1) # A + B <= 1; y = 1 - x
gp = gp + annotate(geom="text", label="A + B > 1\nNot allowed", x=1.02, y=0.07)
gp = gp + geom_abline(slope=1/10, intercept=0) # A <= 10B; y = 0.1x
gp = gp + annotate(geom="text", label="A < 10B\nNot allowed", x=0.35, y=0.05)
gp = gp + geom_hline(yintercept=0.01)
gp = gp + annotate(geom="text", label=expression("B">=0.01), x=1.06, y=0.012)
gp = gp + geom_label(
  data=start_end, mapping=aes(label=Label), alpha=0.7, color="black"
)
gp = gp + xlim(NA, 1.08)
gp = gp + theme_bw()
gp = gp + scale_color_viridis_c()
gp = gp + theme(
  axis.ticks = element_line(colour="black"),
  axis.text = element_text(colour="black"),
  aspect.ratio = 1
)

ggsave("results/toy.sampling.pdf", gp, w=6.5, h=5.5)
