## Figure 2: Nature-style 8-panel composite of Phase 2 (post-fixation) results,
## built entirely from the n = 20 replicated simulation output. Replaces the
## single boxplot with a richer, mechanistically-connected panel set.
library(ggplot2)
library(dplyr)
library(patchwork)

dir.create("../figures", showWarnings = FALSE)

okabe_ito <- c(A = "#0072B2", B = "#009E73", C = "#D55E00", diff = "#CC79A7", acc = "#56B4E9")

theme_nature <- theme_classic(base_size = 9, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    legend.position = "none",
    plot.title = element_blank(),
    plot.margin = margin(4, 8, 4, 4)
  )

p1 <- read.csv("../outputs/phase1_replicated_log.csv")
p1$p_major <- pmax(p1$p_major, 1 - p1$p_major)
p2 <- read.csv("../outputs/phase2_replicated_summary.csv")

gen_fixed_fixed <- sapply(split(p1, p1$rep), function(sub) {
  fg <- sub$generation[sub$p_major >= 0.99]
  if (length(fg) > 0) min(fg) else max(sub$generation)
})
p2$gen_fixed <- gen_fixed_fixed[as.character(p2$rep)]
p2$VA_background_at_fixation <- mapply(function(r, g) {
  p1$VA_background[p1$rep == r & p1$generation == g]
}, p2$rep, p2$gen_fixed)

## ---- A: cumulative gain by strategy (boxplot + jitter) ----
p2_long <- data.frame(
  strategy = rep(c("A", "B", "C"), each = nrow(p2)),
  gain = c(p2$gain_A, p2$gain_B, p2$gain_C)
)
panelA <- ggplot(p2_long, aes(x = strategy, y = gain, fill = strategy, colour = strategy)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.5, linewidth = 0.35) +
  geom_jitter(width = 0.08, size = 0.9, alpha = 0.6) +
  scale_fill_manual(values = unname(okabe_ito[c("A","B","C")])) +
  scale_colour_manual(values = unname(okabe_ito[c("A","B","C")])) +
  labs(x = NULL, y = "Cumulative gain") +
  theme_nature

## ---- B: paired per-replicate trajectories A -> B -> C ----
pspag <- data.frame(rep = rep(p2$rep, 3),
                     strategy = rep(c("A","B","C"), each = nrow(p2)),
                     gain = c(p2$gain_A, p2$gain_B, p2$gain_C))
pspag$strategy <- factor(pspag$strategy, levels = c("A","B","C"))
panelB <- ggplot(pspag, aes(x = strategy, y = gain, group = rep)) +
  geom_line(colour = "grey70", linewidth = 0.3) +
  geom_point(aes(colour = strategy), size = 1) +
  scale_colour_manual(values = unname(okabe_ito[c("A","B","C")])) +
  labs(x = NULL, y = "Gain per replicate (n = 20, paired)") +
  theme_nature

## ---- C: paired difference B - A ----
dAB <- data.frame(d = p2$gain_B - p2$gain_A)
panelC <- ggplot(dAB, aes(x = d)) +
  geom_histogram(bins = 10, fill = unname(okabe_ito["diff"]), colour = "white", linewidth = 0.2, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.35) +
  geom_vline(xintercept = mean(dAB$d), colour = unname(okabe_ito["diff"]), linewidth = 0.6) +
  labs(x = "Gain(B) − Gain(A) per replicate", y = "Count") +
  theme_nature

## ---- D: paired difference B - C ----
dBC <- data.frame(d = p2$gain_B - p2$gain_C)
panelD <- ggplot(dBC, aes(x = d)) +
  geom_histogram(bins = 10, fill = unname(okabe_ito["diff"]), colour = "white", linewidth = 0.2, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.35) +
  geom_vline(xintercept = mean(dBC$d), colour = unname(okabe_ito["diff"]), linewidth = 0.6) +
  labs(x = "Gain(B) − Gain(C) per replicate", y = "Count") +
  theme_nature

## ---- E: genomic-prediction accuracy distribution ----
panelE <- ggplot(p2, aes(x = gs_accuracy)) +
  geom_histogram(bins = 10, fill = unname(okabe_ito["acc"]), colour = "white", linewidth = 0.2, alpha = 0.9) +
  geom_vline(xintercept = mean(p2$gs_accuracy), colour = unname(okabe_ito["acc"]), linewidth = 0.6) +
  labs(x = "GS mean prediction accuracy (r)", y = "Count") +
  theme_nature

## ---- F: generation-to-fixation vs residual background VA (mechanistic link) ----
panelF <- ggplot(p2, aes(x = gen_fixed, y = VA_background_at_fixation)) +
  geom_point(colour = unname(okabe_ito["A"]), size = 1.3, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.5, fill = "grey80") +
  labs(x = "Generation to major-locus fixation", y = expression("Background "*V[A]*" at fixation")) +
  theme_nature

## ---- G: GS accuracy vs realized genomic-selection gain ----
panelG <- ggplot(p2, aes(x = gs_accuracy, y = gain_C)) +
  geom_point(colour = unname(okabe_ito["C"]), size = 1.3, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", linewidth = 0.5, fill = "grey80") +
  labs(x = "GS mean prediction accuracy (r)", y = "Cumulative gain, Arm C") +
  theme_nature

## ---- H: relative advantage of B, with 95% CI from paired t-tests ----
ttAB <- t.test(p2$gain_B, p2$gain_A, paired = TRUE)  ## CI on mean(gain_B - gain_A)
ttBC <- t.test(p2$gain_B, p2$gain_C, paired = TRUE)  ## CI on mean(gain_B - gain_C)
rel <- data.frame(
  comparison = c("B vs A", "B vs C"),
  pct = c(mean(p2$gain_B - p2$gain_A) / mean(p2$gain_A) * 100,
          mean(p2$gain_B - p2$gain_C) / mean(p2$gain_C) * 100),
  lo = c(ttAB$conf.int[1] / mean(p2$gain_A) * 100,
         ttBC$conf.int[1] / mean(p2$gain_C) * 100),
  hi = c(ttAB$conf.int[2] / mean(p2$gain_A) * 100,
         ttBC$conf.int[2] / mean(p2$gain_C) * 100)
)
panelH <- ggplot(rel, aes(x = comparison, y = pct)) +
  geom_col(fill = unname(okabe_ito["B"]), width = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.4) +
  labs(x = NULL, y = "Relative gain advantage of B (%)") +
  theme_nature

fig2 <- (panelA + labs(tag = "A")) + (panelB + labs(tag = "B")) +
        (panelC + labs(tag = "C")) + (panelD + labs(tag = "D")) +
        (panelE + labs(tag = "E")) + (panelF + labs(tag = "F")) +
        (panelG + labs(tag = "G")) + (panelH + labs(tag = "H")) +
  plot_layout(ncol = 2, nrow = 4) &
  theme(plot.tag = element_text(face = "bold", size = 10))

ggsave("../figures/Figure2_three_strategies_comparison.png", fig2, width = 8.5, height = 13, dpi = 320)
ggsave("../figures/Figure2_three_strategies_comparison.pdf", fig2, width = 8.5, height = 13)

cat("8-panel Figure 2 written to simulation/figures/\n")
cat(sprintf("B-A CI (pct pts): [%.2f, %.2f]\n", rel$lo[1], rel$hi[1]))
cat(sprintf("B-C CI (pct pts): [%.2f, %.2f]\n", rel$lo[2], rel$hi[2]))
