## Phase 1: historical fixation of a major-effect quality locus (sh2-like) in a
## simulated sweet-corn population, validating VA = 2pq*a^2 -> 0 as p -> 1.
## Author: N.T. Duc / simulation companion to "Chon giong lay kieu hinh lam trung
## tam trong ky nguyen gen: bai hoc tu ngo thuc pham" (VJAS review).
library(AlphaSimR)

set.seed(20260810)

N_FOUNDER   <- 800
N_CHR       <- 10
SEG_SITES   <- 400
N_BG_QTL_CHR <- 30          # background (polygenic) QTL per chromosome
SEL_INTENSITY <- 0.20        # top 20% selected as parents each generation
N_GEN_PHASE1  <- 18          # ~ decades of commercial fixation, discretised as generations
H2_COMPOSITE  <- 0.35         # heritability of the single composite quality score

## ---- founder haplotypes (realistic coalescent LD via runMacs) ----
founderPop <- runMacs(nInd = N_FOUNDER, nChr = N_CHR, segSites = SEG_SITES,
                       species = "MAIZE", inbred = TRUE)

SP <- SimParam$new(founderPop)
SP$restrSegSites(minQtlPerChr = N_BG_QTL_CHR + 1, minSnpPerChr = 50, overlap = FALSE)

## exactly one major-effect locus, isolated to chromosome 1
major_pos <- rep(0, N_CHR); major_pos[1] <- 1
SP$addTraitA(nQtlPerChr = major_pos, mean = 0, var = 1, name = "major")
## polygenic background: many small-effect loci spread across all chromosomes
SP$addTraitA(nQtlPerChr = rep(N_BG_QTL_CHR, N_CHR), mean = 0, var = 1, name = "background")
SP$addSnpChip(nSnpPerChr = 80)

pop <- newPop(founderPop, simParam = SP)

## composite quality genetic value = major + background (equal a priori weight,
## matches the review's own framing: "chi so chon loc tich hop cam quan")
composite_gv <- function(p) rowSums(gv(p))

logRows <- list()
alleleFreqTraj <- numeric(N_GEN_PHASE1 + 1)
qtlGeno0 <- pullQtlGeno(pop, trait = 1, simParam = SP)
alleleFreqTraj[1] <- mean(qtlGeno0) / 2

record_gen <- function(p, gen) {
  gvmat <- gv(p)
  qtlG <- pullQtlGeno(p, trait = 1, simParam = SP)
  pfreq <- mean(qtlG) / 2
  data.frame(
    generation = gen,
    p_major = pfreq,
    VA_major_observed = var(gvmat[, "major"]),
    VA_background = var(gvmat[, "background"]),
    VA_total = var(rowSums(gvmat)),
    mean_composite = mean(rowSums(gvmat))
  )
}

logRows[[1]] <- record_gen(pop, 0)

for (gen in seq_len(N_GEN_PHASE1)) {
  pheno_composite <- composite_gv(pop) + rnorm(pop@nInd, 0,
                        sqrt(var(composite_gv(pop)) * (1 - H2_COMPOSITE) / H2_COMPOSITE))
  ord <- order(pheno_composite, decreasing = TRUE)
  keepN <- max(20, round(pop@nInd * SEL_INTENSITY))
  parents <- pop[ord[seq_len(keepN)]]
  pop <- randCross(parents, nCrosses = N_FOUNDER, simParam = SP)
  logRows[[gen + 1]] <- record_gen(pop, gen)
}

phase1_log <- do.call(rbind, logRows)

## theoretical VA_major curve from the *actual* major-locus additive effect used
a_major <- SP$traits[[1]]@addEff[1]
phase1_log$VA_major_theory <- with(phase1_log, 2 * p_major * (1 - p_major) * a_major^2)

dir.create("../outputs", showWarnings = FALSE)
write.csv(phase1_log, "../outputs/phase1_fixation_log.csv", row.names = FALSE)
saveRDS(pop, "../outputs/phase1_final_population.rds")
saveRDS(SP, "../outputs/simparam.rds")

cat("Phase 1 complete.\n")
cat("Final p_major:", tail(phase1_log$p_major, 1), "\n")
cat("Final VA_major (observed):", tail(phase1_log$VA_major_observed, 1),
    " | (theory):", tail(phase1_log$VA_major_theory, 1), "\n")
cat("Final VA_background:", tail(phase1_log$VA_background, 1), "\n")
print(phase1_log)
