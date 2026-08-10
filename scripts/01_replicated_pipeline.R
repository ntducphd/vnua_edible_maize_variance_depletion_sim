## Replicated simulation pipeline: N_REP independent replicates, each running
## Phase 1 (historical major-locus fixation) then Phase 2 (three post-fixation
## breeding strategies), for a properly powered statistical comparison.
##
## Note on reproducibility: AlphaSimR's runMacs() calls an external coalescent
## sampler (MaCS) whose own internal RNG is not fully pinned by R's set.seed().
## Per-replicate results therefore vary even at a fixed R seed; this is expected
## and is exactly why replication (rather than a single run) is required for a
## defensible comparison. Aggregate statistics across N_REP independent
## replicates are the reported result, not any single run.
library(AlphaSimR)
library(AGHmatrix)
library(sommer)

N_REP         <- 20
N_FOUNDER     <- 800
N_CHR         <- 10
SEG_SITES     <- 400
N_BG_QTL_CHR  <- 30
SEL_INTENSITY <- 0.20
N_GEN_PHASE1  <- 18
N_GEN_PHASE2  <- 12
H2_COMPOSITE  <- 0.35
H2_SUBTRAIT   <- 0.50
TRAIN_FRAC    <- 0.5

run_one_replicate <- function(rep_id) {
  set.seed(20260810 + rep_id)

  founderPop <- runMacs(nInd = N_FOUNDER, nChr = N_CHR, segSites = SEG_SITES,
                         species = "MAIZE", inbred = TRUE)
  SP <- SimParam$new(founderPop)
  SP$restrSegSites(minQtlPerChr = N_BG_QTL_CHR + 1, minSnpPerChr = 50, overlap = FALSE)
  major_pos <- rep(0, N_CHR); major_pos[1] <- 1
  SP$addTraitA(nQtlPerChr = major_pos, mean = 0, var = 1, name = "major")
  SP$addTraitA(nQtlPerChr = rep(N_BG_QTL_CHR, N_CHR), mean = 0, var = 1, name = "background")
  SP$addSnpChip(nSnpPerChr = 80)

  pop <- newPop(founderPop, simParam = SP)
  composite_gv <- function(p) rowSums(gv(p))

  p1log <- vector("list", N_GEN_PHASE1 + 1)
  qtlG0 <- pullQtlGeno(pop, trait = 1, simParam = SP)
  rec <- function(p, gen) {
    gvmat <- gv(p)
    qtlG <- pullQtlGeno(p, trait = 1, simParam = SP)
    pfreq <- mean(qtlG) / 2
    data.frame(rep = rep_id, generation = gen, p_major = pfreq,
               VA_major = var(gvmat[, "major"]), VA_background = var(gvmat[, "background"]))
  }
  p1log[[1]] <- rec(pop, 0)
  for (g in seq_len(N_GEN_PHASE1)) {
    tc <- composite_gv(pop)
    errSD <- sqrt(var(tc) * (1 - H2_COMPOSITE) / H2_COMPOSITE)
    pheno <- tc + rnorm(pop@nInd, 0, errSD)
    ord <- order(pheno, decreasing = TRUE)
    keepN <- max(20, round(pop@nInd * SEL_INTENSITY))
    parents <- pop[ord[seq_len(keepN)]]
    pop <- randCross(parents, nCrosses = N_FOUNDER, simParam = SP)
    p1log[[g + 1]] <- rec(pop, g)
  }
  phase1_df <- do.call(rbind, p1log)
  a_major <- SP$traits[[1]]@addEff[1]
  phase1_df$VA_major_theory <- with(phase1_df, 2 * p_major * (1 - p_major) * a_major^2)
  fixed_gens <- phase1_df$generation[phase1_df$p_major >= 0.99]
  ## fall back to the final Phase-1 generation if the locus has not reached
  ## p >= 0.99 within N_GEN_PHASE1 (a real, occasional stochastic outcome,
  ## not an error) -- Phase 2 still starts from that final population either way
  gen_fixed <- if (length(fixed_gens) > 0) min(fixed_gens) else N_GEN_PHASE1

  ## ---- Phase 2: three arms from this replicate's fixed population ----
  pop0 <- pop
  bgEff <- SP$traits[[2]]@addEff
  grp <- sample(rep(1:4, length.out = length(bgEff)))
  true_composite <- function(p) rowSums(gv(p))
  sub_scores <- function(p) {
    qtlG <- pullQtlGeno(p, trait = 2, simParam = SP)
    sapply(1:4, function(k) as.numeric(qtlG[, grp == k, drop = FALSE] %*% bgEff[grp == k]))
  }

  run_arm_A <- function(pop, nGen) {
    m0 <- mean(true_composite(pop))
    for (g in seq_len(nGen)) {
      tc <- true_composite(pop)
      errSD <- sqrt(var(tc) * (1 - H2_COMPOSITE) / H2_COMPOSITE)
      pheno <- tc + rnorm(pop@nInd, 0, errSD)
      ord <- order(pheno, decreasing = TRUE)
      keepN <- max(20, round(pop@nInd * SEL_INTENSITY))
      pop <- randCross(pop[ord[seq_len(keepN)]], nCrosses = pop0@nInd, simParam = SP)
    }
    mean(true_composite(pop)) - m0
  }

  run_arm_B <- function(pop, nGen) {
    m0 <- mean(true_composite(pop))
    for (g in seq_len(nGen)) {
      sc <- sub_scores(pop)
      errSD_k <- sqrt(apply(sc, 2, var) * (1 - H2_SUBTRAIT) / H2_SUBTRAIT)
      phenoK <- sapply(1:4, function(k) sc[, k] + rnorm(pop@nInd, 0, errSD_k[k]))
      VPk <- apply(phenoK, 2, var); VAk <- apply(sc, 2, var)
      index <- as.numeric(phenoK %*% (VAk / VPk))
      ord <- order(index, decreasing = TRUE)
      keepN <- max(20, round(pop@nInd * SEL_INTENSITY))
      pop <- randCross(pop[ord[seq_len(keepN)]], nCrosses = pop0@nInd, simParam = SP)
    }
    mean(true_composite(pop)) - m0
  }

  run_arm_C <- function(pop, nGen) {
    m0 <- mean(true_composite(pop))
    accs <- numeric(nGen)
    for (g in seq_len(nGen)) {
      tc <- true_composite(pop)
      errSD <- sqrt(var(tc) * (1 - H2_COMPOSITE) / H2_COMPOSITE)
      pheno <- tc + rnorm(pop@nInd, 0, errSD)
      M <- pullSnpGeno(pop, simParam = SP)
      G <- Gmatrix(M, method = "VanRaden", ploidy = 2) + diag(1e-4, pop@nInd)
      trainIdx <- sample(seq_len(pop@nInd), round(pop@nInd * TRAIN_FRAC))
      y <- rep(NA_real_, pop@nInd); y[trainIdx] <- pheno[trainIdx]
      dat <- data.frame(id = as.character(seq_len(pop@nInd)), y = y)
      rownames(G) <- colnames(G) <- dat$id
      fit <- mmer(y ~ 1, random = ~ vsr(id, Gu = G), rcov = ~ units, data = dat, verbose = FALSE)
      gebv <- fit$U$`u:id`$y[dat$id]
      accs[g] <- suppressWarnings(cor(gebv[-trainIdx], tc[-trainIdx]))
      ord <- order(gebv, decreasing = TRUE)
      keepN <- max(20, round(pop@nInd * SEL_INTENSITY))
      pop <- randCross(pop[ord[seq_len(keepN)]], nCrosses = pop0@nInd, simParam = SP)
    }
    c(gain = mean(true_composite(pop)) - m0, acc = mean(accs, na.rm = TRUE))
  }

  gainA <- run_arm_A(pop0, N_GEN_PHASE2)
  gainB <- run_arm_B(pop0, N_GEN_PHASE2)
  resC  <- run_arm_C(pop0, N_GEN_PHASE2)

  phase2_df <- data.frame(rep = rep_id,
                           gain_A = gainA, gain_B = gainB,
                           gain_C = resC["gain"], gs_accuracy = resC["acc"],
                           gen_fixed = gen_fixed,
                           VA_background_at_fixation = phase1_df$VA_background[phase1_df$generation == gen_fixed],
                           VA_background_final = tail(phase1_df$VA_background, 1))

  list(phase1 = phase1_df, phase2 = phase2_df)
}

cat("Running", N_REP, "independent replicates (Phase 1 + Phase 2 each)...\n")
dir.create("../outputs", showWarnings = FALSE)
p1_path <- "../outputs/phase1_replicated_log.csv"
p2_path <- "../outputs/phase2_replicated_summary.csv"
if (file.exists(p1_path)) file.remove(p1_path)
if (file.exists(p2_path)) file.remove(p2_path)

all_p1 <- vector("list", N_REP)
all_p2 <- vector("list", N_REP)
for (r in seq_len(N_REP)) {
  t0 <- Sys.time()
  res <- tryCatch(run_one_replicate(r), error = function(e) {
    cat(sprintf("Replicate %d FAILED: %s -- retrying with a new seed offset\n", r, conditionMessage(e)))
    NULL
  })
  if (is.null(res)) {
    set.seed(90000 + r)
    res <- run_one_replicate(r)
  }
  all_p1[[r]] <- res$phase1
  all_p2[[r]] <- res$phase2
  ## incremental save so a later crash does not lose completed replicates
  write.table(res$phase1, p1_path, sep = ",", row.names = FALSE,
              col.names = !file.exists(p1_path), append = file.exists(p1_path))
  write.table(res$phase2, p2_path, sep = ",", row.names = FALSE,
              col.names = !file.exists(p2_path), append = file.exists(p2_path))
  cat(sprintf("Replicate %d/%d done in %.1fs | fixed at gen %d | gains A=%.3f B=%.3f C=%.3f | GSacc=%.3f\n",
              r, N_REP, as.numeric(Sys.time() - t0, units = "secs"),
              res$phase2$gen_fixed, res$phase2$gain_A, res$phase2$gain_B,
              res$phase2$gain_C, res$phase2$gs_accuracy))
}

phase1_all <- do.call(rbind, all_p1)
phase2_all <- do.call(rbind, all_p2)

cat("\n=== Summary across", N_REP, "replicates ===\n")
cat(sprintf("Generation to fixation: mean=%.2f sd=%.2f\n",
            mean(phase2_all$gen_fixed), sd(phase2_all$gen_fixed)))
cat(sprintf("Cumulative gain, Arm A (phenotypic): mean=%.3f sd=%.3f\n",
            mean(phase2_all$gain_A), sd(phase2_all$gain_A)))
cat(sprintf("Cumulative gain, Arm B (multi-trait index): mean=%.3f sd=%.3f\n",
            mean(phase2_all$gain_B), sd(phase2_all$gain_B)))
cat(sprintf("Cumulative gain, Arm C (genomic selection): mean=%.3f sd=%.3f\n",
            mean(phase2_all$gain_C), sd(phase2_all$gain_C)))
cat(sprintf("GS mean prediction accuracy: mean=%.3f sd=%.3f\n",
            mean(phase2_all$gs_accuracy), sd(phase2_all$gs_accuracy)))

cat("\nPaired t-tests (two-sided):\n")
print(t.test(phase2_all$gain_B, phase2_all$gain_A, paired = TRUE))
print(t.test(phase2_all$gain_B, phase2_all$gain_C, paired = TRUE))
print(t.test(phase2_all$gain_A, phase2_all$gain_C, paired = TRUE))
