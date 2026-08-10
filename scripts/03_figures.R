## Figures from real simulation output (no invented numbers).
## Fig S1: VA_major -> 0 as p -> 1 (Phase 1, validates VA = 2pq*a^2)
## Fig S2: cumulative genetic gain, 3 post-fixation strategies (Phase 2)
library(ggplot2)

dir.create("../figures", showWarnings = FALSE)

## ---- Fig S1 ----
p1log <- read.csv("../outputs/phase1_fixation_log.csv")

fig_s1a <- ggplot(p1log, aes(x = generation)) +
  geom_line(aes(y = p_major), linewidth = 1, color = "#D55E00") +
  geom_point(aes(y = p_major), color = "#D55E00", size = 2) +
  labs(x = "The he chon loc", y = "Tan so alen thuan loi tai locus lon (p)") +
  scale_x_continuous(breaks = seq(0, max(p1log$generation), 2)) +
  ylim(0, 1) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

fig_s1b <- ggplot(p1log, aes(x = generation)) +
  geom_line(aes(y = VA_major_observed, color = "Gen lon (VA quan sat)"), linewidth = 1) +
  geom_point(aes(y = VA_major_observed, color = "Gen lon (VA quan sat)"), size = 2) +
  geom_line(aes(y = VA_major_theory, color = "Gen lon (ly thuyet 2pq.a2)"),
            linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = VA_background, color = "Nen da gen"), linewidth = 1) +
  geom_point(aes(y = VA_background, color = "Nen da gen"), size = 2) +
  scale_color_manual(values = c("Gen lon (VA quan sat)" = "#D55E00",
                                  "Gen lon (ly thuyet 2pq.a2)" = "#D55E00",
                                  "Nen da gen" = "#0072B2"),
                      name = NULL) +
  labs(x = "The he chon loc", y = "Phuong sai di truyen cong gop (VA)") +
  scale_x_continuous(breaks = seq(0, max(p1log$generation), 2)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("../figures/figS1a_allele_freq.png", fig_s1a, width = 6.5, height = 4.5, dpi = 300)
ggsave("../figures/figS1b_VA_depletion.png", fig_s1b, width = 6.5, height = 5, dpi = 300)

## ---- Fig S2 ----
p2log <- read.csv("../outputs/phase2_three_strategies_log.csv")
arm_labels <- c(
  A_phenotypic_mass = "Chon loc kieu hinh don tinh trang",
  B_multitrait_index = "Chon loc chi so da tinh trang",
  C_genomic_selection = "Chon loc bo gen (GBLUP)"
)
p2log$arm_label <- arm_labels[p2log$arm]

fig_s2 <- ggplot(p2log, aes(x = generation, y = mean_true, color = arm_label,
                              shape = arm_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("Chon loc kieu hinh don tinh trang" = "#0072B2",
                                  "Chon loc chi so da tinh trang" = "#009E73",
                                  "Chon loc bo gen (GBLUP)" = "#D55E00"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  labs(x = "The he sau khi gen lon da co dinh",
       y = "Gia tri di truyen tong hop trung binh (chi so chat luong)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave("../figures/figS2_three_strategies_gain.png", fig_s2, width = 7, height = 5, dpi = 300)

cat("Figures written to simulation/figures/\n")
