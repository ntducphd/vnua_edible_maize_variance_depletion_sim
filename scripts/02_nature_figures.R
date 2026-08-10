## Nature-style figures from the replicated simulation (N = 20 independent
## replicates). Okabe-Ito colourblind-safe palette, panel labels, no baked-in
## titles, ribbons/error bars derived from real replicate variance.
library(ggplot2)
library(dplyr)

dir.create("../figures", showWarnings = FALSE)

okabe_ito <- c(major = "#D55E00", background = "#0072B2",
               A = "#0072B2", B = "#009E73", C = "#D55E00")

theme_nature <- theme_classic(base_size = 17, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks = element_line(linewidth = 0.5, colour = "black"),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 15, hjust = 0),
    plot.title = element_blank()
  )

## ---------------- Figure 1: Phase 1 (fixation) ----------------
p1 <- read.csv("../outputs/phase1_replicated_log.csv")
## see note in 03_tables.R: report allele-agnostic distance to fixation
p1$p_major <- pmax(p1$p_major, 1 - p1$p_major)

p1_summary <- p1 %>%
  group_by(generation) %>%
  summarise(
    p_mean = mean(p_major), p_sd = sd(p_major),
    VAmaj_mean = mean(VA_major), VAmaj_sd = sd(VA_major),
    VAbg_mean = mean(VA_background), VAbg_sd = sd(VA_background),
    VAtheory_mean = mean(VA_major_theory),
    .groups = "drop"
  )

figA <- ggplot(p1_summary, aes(x = generation)) +
  geom_ribbon(aes(ymin = pmax(0, p_mean - p_sd), ymax = pmin(1, p_mean + p_sd)),
              fill = okabe_ito["major"], alpha = 0.2) +
  geom_line(aes(y = p_mean), colour = okabe_ito["major"], linewidth = 0.9) +
  geom_point(aes(y = p_mean), colour = okabe_ito["major"], size = 1.6) +
  scale_x_continuous(breaks = seq(0, max(p1_summary$generation), 3)) +
  labs(x = "Selection generation", y = "Favourable allele frequency (p)") +
  theme_nature

figB <- ggplot(p1_summary, aes(x = generation)) +
  geom_ribbon(aes(ymin = pmax(0, VAmaj_mean - VAmaj_sd), ymax = VAmaj_mean + VAmaj_sd,
                  fill = "Major locus"), alpha = 0.18) +
  geom_ribbon(aes(ymin = pmax(0, VAbg_mean - VAbg_sd), ymax = VAbg_mean + VAbg_sd,
                  fill = "Polygenic background"), alpha = 0.18) +
  geom_line(aes(y = VAmaj_mean, colour = "Major locus (observed)"), linewidth = 0.9) +
  geom_line(aes(y = VAtheory_mean, colour = "Major locus (theory, 2pq·a²)"),
            linetype = "dashed", linewidth = 0.7) +
  geom_line(aes(y = VAbg_mean, colour = "Polygenic background"), linewidth = 0.9) +
  scale_colour_manual(values = c("Major locus (observed)" = unname(okabe_ito["major"]),
                                   "Major locus (theory, 2pq·a²)" = unname(okabe_ito["major"]),
                                   "Polygenic background" = unname(okabe_ito["background"]))) +
  scale_fill_manual(values = c("Major locus" = unname(okabe_ito["major"]),
                                 "Polygenic background" = unname(okabe_ito["background"])), guide = "none") +
  scale_x_continuous(breaks = seq(0, max(p1_summary$generation), 3)) +
  labs(x = "Selection generation", y = expression("Additive genetic variance ("*V[A]*")")) +
  theme_nature + theme(legend.text = element_text(size = 14))

library(patchwork)
fig1 <- (figA + labs(tag = "A")) + (figB + labs(tag = "B")) +
  plot_layout(ncol = 2) &
  theme(plot.tag = element_text(face = "bold", size = 20))

ggsave("../figures/Figure1_major_locus_fixation.png", fig1, width = 12, height = 5.6, dpi = 320)
ggsave("../figures/Figure1_major_locus_fixation.pdf", fig1, width = 12, height = 5.6)

## ---------------- Figure 2: Phase 2 (three strategies) ----------------
p2 <- read.csv("../outputs/phase2_replicated_summary.csv")
p2_long <- data.frame(
  strategy = rep(c("A. Phenotypic\nmass selection", "B. Multi-trait\nindex selection",
                    "C. Genomic\nselection (GBLUP)"), each = nrow(p2)),
  gain = c(p2$gain_A, p2$gain_B, p2$gain_C)
)
p2_long$strategy <- factor(p2_long$strategy, levels = unique(p2_long$strategy))

fig2 <- ggplot(p2_long, aes(x = strategy, y = gain, fill = strategy)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.08, size = 1.4, alpha = 0.6, aes(colour = strategy)) +
  scale_fill_manual(values = unname(c(okabe_ito["A"], okabe_ito["B"], okabe_ito["C"])), guide = "none") +
  scale_colour_manual(values = unname(c(okabe_ito["A"], okabe_ito["B"], okabe_ito["C"])), guide = "none") +
  labs(x = NULL, y = "Cumulative genetic gain, 12 post-fixation generations\n(n = 20 independent replicates)") +
  theme_nature + theme(legend.position = "none")

ggsave("../figures/Figure2_three_strategies_comparison.png", fig2, width = 6, height = 5.2, dpi = 300)
ggsave("../figures/Figure2_three_strategies_comparison.pdf", fig2, width = 6, height = 5.2)

cat("Nature-style figures written to simulation/figures/\n")
