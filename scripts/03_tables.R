## Build the three manuscript tables from real replicated output.
library(dplyr)

dir.create("../outputs/tables", showWarnings = FALSE)

## ---------------- Table 1: simulation parameters (static) ----------------
table1 <- data.frame(
  Parameter = c("Founder population size", "Chromosomes", "Segregating sites per chromosome",
                "Major-effect QTL", "Background QTL (total)", "Background QTL per chromosome",
                "SNP markers (genomic prediction panel)", "Selection intensity (top fraction retained)",
                "Phase 1 generations (historical fixation)", "Phase 2 generations (post-fixation)",
                "Heritability, composite phenotype (Arms A, C)", "Heritability, sub-trait phenotype (Arm B)",
                "Genomic-selection training fraction per generation", "Independent replicates"),
  Value = c("800", "10", "400", "1 (isolated to chromosome 1)", "300", "30",
            "800", "20%", "18", "12", "0.35", "0.50", "50%", "20")
)
write.csv(table1, "../outputs/tables/Table1_parameters.csv", row.names = FALSE)

## ---------------- Table 2: Phase 1 fixation summary ----------------
p1 <- read.csv("../outputs/phase1_replicated_log.csv")
p2s <- read.csv("../outputs/phase2_replicated_summary.csv")
## AlphaSimR assigns the major-locus effect sign to the reference allele
## arbitrarily per replicate, so raw p_major drifts toward 1 in some
## replicates and toward 0 in others (both represent fixation of whichever
## allele is under positive selection). Report the allele-agnostic distance
## to fixation instead, which is what VA = 2pq*a^2 actually depends on.
p1$p_major <- pmax(p1$p_major, 1 - p1$p_major)

## gen_fixed and VA_background_at_fixation stored in phase2_replicated_summary.csv
## were computed inside the simulation run using the *raw* (asymmetric) p_major,
## so they inherit the same bug for replicates that fixed toward 0 -- recompute
## both properly here from the corrected per-replicate trajectories.
gen_fixed_fixed <- sapply(split(p1, p1$rep), function(sub) {
  fg <- sub$generation[sub$p_major >= 0.99]
  if (length(fg) > 0) min(fg) else max(sub$generation)
})
p2s$gen_fixed <- gen_fixed_fixed[as.character(p2s$rep)]
p2s$VA_background_at_fixation <- mapply(function(r, g) {
  p1$VA_background[p1$rep == r & p1$generation == g]
}, p2s$rep, p2s$gen_fixed)

checkpoints <- c(0, 3, 6, 9, 12, 15, 18)
t2rows <- lapply(checkpoints, function(g) {
  sub <- p1[p1$generation == g, ]
  data.frame(
    Generation = g,
    `Allele freq (p), mean +/- SD` = sprintf("%.3f +/- %.3f", mean(sub$p_major), sd(sub$p_major)),
    `VA major locus, mean +/- SD` = sprintf("%.3f +/- %.3f", mean(sub$VA_major), sd(sub$VA_major)),
    `VA background, mean +/- SD` = sprintf("%.3f +/- %.3f", mean(sub$VA_background), sd(sub$VA_background)),
    check.names = FALSE
  )
})
table2 <- do.call(rbind, t2rows)
write.csv(table2, "../outputs/tables/Table2_phase1_fixation.csv", row.names = FALSE)

## ---------------- Table 3: Phase 2 strategy comparison ----------------
ttestAB <- t.test(p2s$gain_B, p2s$gain_A, paired = TRUE)
ttestBC <- t.test(p2s$gain_B, p2s$gain_C, paired = TRUE)
ttestAC <- t.test(p2s$gain_A, p2s$gain_C, paired = TRUE)

fmt_p <- function(p) if (p < 0.001) "< 0.001" else sprintf("%.3f", p)

table3 <- data.frame(
  Strategy = c("A. Phenotypic mass selection", "B. Multi-trait index selection",
               "C. Genomic selection (GBLUP)"),
  `Cumulative gain, mean +/- SD` = c(
    sprintf("%.3f +/- %.3f", mean(p2s$gain_A), sd(p2s$gain_A)),
    sprintf("%.3f +/- %.3f", mean(p2s$gain_B), sd(p2s$gain_B)),
    sprintf("%.3f +/- %.3f", mean(p2s$gain_C), sd(p2s$gain_C))
  ),
  `Relative to Arm A` = c("1.00x",
                            sprintf("%.2fx", mean(p2s$gain_B) / mean(p2s$gain_A)),
                            sprintf("%.2fx", mean(p2s$gain_C) / mean(p2s$gain_A))),
  check.names = FALSE
)
write.csv(table3, "../outputs/tables/Table3_strategy_comparison.csv", row.names = FALSE)

table3_stats <- data.frame(
  Comparison = c("B vs A (multi-trait index vs phenotypic mass)",
                 "B vs C (multi-trait index vs genomic selection)",
                 "A vs C (phenotypic mass vs genomic selection)"),
  `Mean difference` = c(mean(p2s$gain_B - p2s$gain_A), mean(p2s$gain_B - p2s$gain_C),
                          mean(p2s$gain_A - p2s$gain_C)),
  `t statistic` = c(ttestAB$statistic, ttestBC$statistic, ttestAC$statistic),
  df = c(ttestAB$parameter, ttestBC$parameter, ttestAC$parameter),
  `p value` = c(fmt_p(ttestAB$p.value), fmt_p(ttestBC$p.value), fmt_p(ttestAC$p.value)),
  check.names = FALSE
)
write.csv(table3_stats, "../outputs/tables/Table3b_paired_ttests.csv", row.names = FALSE)

cat("Tables written to simulation/outputs/tables/\n")
print(table1); cat("\n"); print(table2); cat("\n"); print(table3); cat("\n"); print(table3_stats)

## ---------------- values for MANUSCRIPT.md placeholders ----------------
cat("\n\n=== MANUSCRIPT.md placeholder values ===\n")
cat("GEN_FIXED_MEAN:", sprintf("%.1f", mean(p2s$gen_fixed)), "\n")
cat("GEN_FIXED_SD:", sprintf("%.1f", sd(p2s$gen_fixed)), "\n")
cat("VA_BG_FIXED_MEAN:", sprintf("%.3f", mean(p2s$VA_background_at_fixation)), "\n")
cat("VA_BG_FINAL_MEAN:", sprintf("%.3f", mean(p2s$VA_background_final)), "\n")
cat("VA_MAJOR_INIT:", sprintf("%.3f", mean(p1$VA_major[p1$generation == 0])), "\n")
cat("VA_MAJOR_FINAL:", sprintf("%.4f", mean(p1$VA_major[p1$generation == max(p1$generation)])), "\n")
cat("PCT_B_OVER_A:", sprintf("%.1f", (mean(p2s$gain_B)/mean(p2s$gain_A) - 1) * 100), "\n")
cat("PCT_B_OVER_C:", sprintf("%.1f", (mean(p2s$gain_B)/mean(p2s$gain_C) - 1) * 100), "\n")
cat("ACC_MEAN:", sprintf("%.3f", mean(p2s$gs_accuracy)), "\n")
cat("ACC_SD:", sprintf("%.3f", sd(p2s$gs_accuracy)), "\n")
cat("GAIN_A_MEAN:", sprintf("%.3f", mean(p2s$gain_A)), " GAIN_A_SD:", sprintf("%.3f", sd(p2s$gain_A)), "\n")
cat("GAIN_B_MEAN:", sprintf("%.3f", mean(p2s$gain_B)), " GAIN_B_SD:", sprintf("%.3f", sd(p2s$gain_B)), "\n")
cat("GAIN_C_MEAN:", sprintf("%.3f", mean(p2s$gain_C)), " GAIN_C_SD:", sprintf("%.3f", sd(p2s$gain_C)), "\n")
cat("TTEST_AB: t(", ttestAB$parameter, ")=", sprintf("%.2f", ttestAB$statistic), ", p", fmt_p(ttestAB$p.value), "\n")
cat("TTEST_BC: t(", ttestBC$parameter, ")=", sprintf("%.2f", ttestBC$statistic), ", p", fmt_p(ttestBC$p.value), "\n")
