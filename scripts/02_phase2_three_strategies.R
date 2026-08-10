## Phase 2: three post-fixation breeding strategies competing for the residual
## polygenic (background) variance identified in Phase 1.
##   Arm A - classical single-score phenotypic mass selection
##   Arm B - multi-trait index selection on 4 decomposed quality sub-traits
##           ("vi ngot", "do mem vo", "huong thom", "kha nang giu chat luong")
##   Arm C - genomic selection (GBLUP, VanRaden G-matrix via AGHmatrix + sommer)
library(AlphaSimR)
library(AGHmatrix)
library(sommer)

set.seed(20260811)

pop0 <- readRDS("../outputs/phase1_final_population.rds")
SP   <- readRDS("../outputs/simparam.rds")

N_GEN_PHASE2  <- 12
SEL_INTENSITY <- 0.20
H2_COMPOSITE  <- 0.35    # Arm A: single holistic score
H2_SUBTRAIT   <- 0.50    # Arm B: narrower sub-traits, less measurement noise each
TRAIN_FRAC    <- 0.5     # Arm C: fraction phenotyped as GS training population

bgEff <- SP$traits[[2]]@addEff
nBgQtl <- length(bgEff)
set.seed(1)
grp <- sample(rep(1:4, length.out = nBgQtl))   # 4 disjoint sub-trait groups
subtrait_names <- c("vi_ngot", "do_mem_vo", "huong_thom", "giu_chat_luong")

true_composite <- function(p) {
  gvmat <- gv(p)
  rowSums(gvmat)   # major (now ~0 variance) + background = true breeding objective
}

sub_scores <- function(p) {
  qtlG <- pullQtlGeno(p, trait = 2, simParam = SP)
  sapply(1:4, function(k) as.numeric(qtlG[, grp == k, drop = FALSE] %*% bgEff[grp == k]))
}

run_one_gen_gain <- function(p) {
  tc <- true_composite(p)
  c(mean = mean(tc), var = var(tc))
}

## ---------- Arm A: phenotypic mass selection ----------
run_arm_A <- function(pop, nGen) {
  log <- data.frame(generation = 0, mean_true = mean(true_composite(pop)),
                     var_true = var(true_composite(pop)))
  for (g in seq_len(nGen)) {
    tc <- true_composite(pop)
    errSD <- sqrt(var(tc) * (1 - H2_COMPOSITE) / H2_COMPOSITE)
    pheno <- tc + rnorm(pop@nInd, 0, errSD)
    ord <- order(pheno, decreasing = TRUE)
    keepN <- max(20, round(pop@nInd * SEL_INTENSITY))
    parents <- pop[ord[seq_len(keepN)]]
    pop <- randCross(parents, nCrosses = pop0@nInd, simParam = SP)
    log <- rbind(log, data.frame(generation = g, mean_true = mean(true_composite(pop)),
                                  var_true = var(true_composite(pop))))
  }
  list(pop = pop, log = log)
}

## ---------- Arm B: multi-trait index selection ----------
run_arm_B <- function(pop, nGen) {
  log <- data.frame(generation = 0, mean_true = mean(true_composite(pop)),
                     var_true = var(true_composite(pop)))
  for (g in seq_len(nGen)) {
    sc <- sub_scores(pop)                       # true genetic values per sub-trait
    errSD_k <- sqrt(apply(sc, 2, var) * (1 - H2_SUBTRAIT) / H2_SUBTRAIT)
    phenoK <- sapply(1:4, function(k) sc[, k] + rnorm(pop@nInd, 0, errSD_k[k]))
    ## optimal index weight per sub-trait = VA_k / VP_k (regression of true value on
    ## each phenotyped component; sub-traits are genetically independent by
    ## construction so weights do not need to account for genetic covariances)
    VPk <- apply(phenoK, 2, var)
    VAk <- apply(sc, 2, var)
    wk  <- VAk / VPk
    index <- as.numeric(phenoK %*% wk)
    ord <- order(index, decreasing = TRUE)
    keepN <- max(20, round(pop@nInd * SEL_INTENSITY))
    parents <- pop[ord[seq_len(keepN)]]
    pop <- randCross(parents, nCrosses = pop0@nInd, simParam = SP)
    log <- rbind(log, data.frame(generation = g, mean_true = mean(true_composite(pop)),
                                  var_true = var(true_composite(pop))))
  }
  list(pop = pop, log = log)
}

## ---------- Arm C: genomic selection (GBLUP) ----------
run_arm_C <- function(pop, nGen) {
  log <- data.frame(generation = 0, mean_true = mean(true_composite(pop)),
                     var_true = var(true_composite(pop)), gs_accuracy = NA)
  for (g in seq_len(nGen)) {
    tc <- true_composite(pop)
    errSD <- sqrt(var(tc) * (1 - H2_COMPOSITE) / H2_COMPOSITE)
    pheno <- tc + rnorm(pop@nInd, 0, errSD)

    M <- pullSnpGeno(pop, simParam = SP)               # n x m, coded 0/1/2
    G <- Gmatrix(M, method = "VanRaden", ploidy = 2)
    G <- G + diag(1e-4, nrow(G))                        # numerical stability

    trainIdx <- sample(seq_len(pop@nInd), round(pop@nInd * TRAIN_FRAC))
    y <- rep(NA_real_, pop@nInd)
    y[trainIdx] <- pheno[trainIdx]

    dat <- data.frame(id = as.character(seq_len(pop@nInd)), y = y)
    rownames(G) <- colnames(G) <- dat$id
    fit <- mmer(y ~ 1, random = ~ vsr(id, Gu = G), rcov = ~ units,
                data = dat, verbose = FALSE)
    gebv <- fit$U$`u:id`$y[dat$id]

    acc <- suppressWarnings(cor(gebv[-trainIdx], tc[-trainIdx]))

    ord <- order(gebv, decreasing = TRUE)
    keepN <- max(20, round(pop@nInd * SEL_INTENSITY))
    parents <- pop[ord[seq_len(keepN)]]
    pop <- randCross(parents, nCrosses = pop0@nInd, simParam = SP)
    log <- rbind(log, data.frame(generation = g, mean_true = mean(true_composite(pop)),
                                  var_true = var(true_composite(pop)), gs_accuracy = acc))
  }
  list(pop = pop, log = log)
}

cat("Starting Arm A (phenotypic mass selection)...\n")
resA <- run_arm_A(pop0, N_GEN_PHASE2)
cat("Starting Arm B (multi-trait index selection)...\n")
resB <- run_arm_B(pop0, N_GEN_PHASE2)
cat("Starting Arm C (genomic selection)...\n")
resC <- run_arm_C(pop0, N_GEN_PHASE2)

resA$log$arm <- "A_phenotypic_mass"
resB$log$arm <- "B_multitrait_index"
resC$log$arm <- "C_genomic_selection"
resC$log$gs_accuracy <- resC$log$gs_accuracy

allLog <- rbind(resA$log[, c("generation","mean_true","var_true","arm")],
                resB$log[, c("generation","mean_true","var_true","arm")],
                resC$log[, c("generation","mean_true","var_true","arm")])

write.csv(allLog, "../outputs/phase2_three_strategies_log.csv", row.names = FALSE)
write.csv(resC$log, "../outputs/phase2_armC_gs_accuracy.csv", row.names = FALSE)

cat("\n=== Cumulative genetic gain after", N_GEN_PHASE2, "cycles ===\n")
for (nm in c("A_phenotypic_mass", "B_multitrait_index", "C_genomic_selection")) {
  sub <- allLog[allLog$arm == nm, ]
  gain <- tail(sub$mean_true, 1) - sub$mean_true[1]
  cat(sprintf("%s: start=%.3f end=%.3f gain=%.3f\n", nm, sub$mean_true[1],
              tail(sub$mean_true,1), gain))
}
cat("\nArm C mean GS accuracy across cycles:", mean(resC$log$gs_accuracy, na.rm = TRUE), "\n")
