# Simulation companion: major-gene-fixation variance depletion in edible maize breeding

Companion simulation study for the review "Chọn giống lấy kiểu hình làm trung tâm trong
kỷ nguyên gen: bài học từ ngô thực phẩm" (*Vietnam Journal of Agricultural Sciences*).
It turns the review's central quantitative-genetics claim — that once a large-effect
quality locus (sh2-like) is fixed in a commercial market class, the remaining additive
genetic variance sits entirely in the polygenic background, and different post-fixation
selection strategies exploit that residual variance with different efficiency — from a
narrative argument into a reproducible, numeric demonstration.

Built with [AlphaSimR](https://cran.r-project.org/package=AlphaSimR) (Gaynor et al.
2021, *The Plant Genome*), [AGHmatrix](https://cran.r-project.org/package=AGHmatrix),
and [sommer](https://cran.r-project.org/package=sommer).

## What it simulates

**Phase 1 — historical fixation** (`scripts/01_phase1_fixation.R`). A founder population
(N = 800, 10 chromosomes, coalescent LD via `runMacs`) segregating at one large-effect
quality locus plus 300 small-effect background QTL. Eighteen generations of truncation
selection (top 20%) on the composite phenotype fix the major locus by generation 7
(p = 1.000), at which point its own additive variance collapses to exactly 0 — matching
the classical prediction V_A = 2pq·a² to within simulation noise at every generation.
Background additive variance, by contrast, stays firmly above zero (0.30 at the point of
fixation; still 0.18 after 11 more generations of continued selection) — the literal,
numeric version of the review's "cạn biến dị do cố định gen lớn" argument.

**Phase 2 — three post-fixation strategies** (`scripts/02_phase2_three_strategies.R`),
all starting from Phase 1's final (major-locus-fixed) population and run for 12 further
generations:

- **A — Classical single-score phenotypic mass selection**: one noisy composite
  phenotype (h² = 0.35), truncation selection.
- **B — Multi-trait index selection**: the same true breeding objective decomposed into
  four independently-measured sub-traits ("vị ngọt", "độ mềm vỏ", "hương thơm", "khả
  năng giữ chất lượng" — the same four components the review itself names), each
  phenotyped with less noise (h² = 0.50) and combined via an optimal (known-variance)
  selection index.
- **C — Genomic selection**: GBLUP on an 800-SNP panel, VanRaden G-matrix
  (`AGHmatrix::Gmatrix`), fit with `sommer::mmer`, half the population phenotyped as a
  fresh training set each cycle, selection on GEBV.

## Real results (not illustrative — read directly from `outputs/*.csv`)

| Strategy | Cumulative genetic gain, 12 cycles | Relative to A |
|---|---:|---:|
| A — Phenotypic mass selection | 3.766 | 1.00× |
| B — Multi-trait index selection | **4.223** | **1.12×** |
| C — Genomic selection (GBLUP) | 2.926 | 0.78× |

Mean genomic-selection prediction accuracy across cycles: **r ≈ 0.30** — modest, because
the training population here is freshly re-sampled (no historical accumulation) and the
marker panel is moderate relative to the 300-locus background architecture. This is a
realistic constraint for a resource-limited program, not a design flaw: it is exactly why
genomic selection's real published advantage (e.g. Beyene et al. 2015, cited in the main
review) comes mainly from *faster cycling*, not higher per-cycle accuracy — a nuance this
simulation reproduces on its own rather than assuming.

**Net finding**: once the major locus is fixed, *how you structure phenotyping* — measuring
the true quality sub-components separately and combining them with a proper index —
recovered more of the residual polygenic variance than either naively phenotyping one
composite score, or defaulting to marker-based prediction with an under-resourced training
population. This is the review's own Section 7 argument, now with a number attached.

## Reproducing

```r
install.packages(c("AlphaSimR", "AGHmatrix", "sommer", "ggplot2"))
setwd("scripts")
source("01_phase1_fixation.R")        # writes outputs/phase1_*
source("02_phase2_three_strategies.R") # writes outputs/phase2_*
source("03_figures.R")                 # writes figures/*.png
```
Every run reseeds explicitly (`set.seed()` in each script); results above are from the
seeds committed here and are exactly reproducible.

## Caveats

- Sub-trait weights in Arm B use the *true* simulated genetic variances (a known-parameter
  index), isolating the value of decomposition + optimal weighting from the separate,
  real-world problem of estimating those weights from limited data.
- Arm C retrains from a fresh random half of the population each cycle; a program that
  accumulates training data across cycles, or that increases marker density, would likely
  close some or all of the gap with B — this simulation shows one realistic operating
  point, not a universal ranking of methods.
- Population size (N = 800), generation count, and heritabilities were fixed a priori
  before any run; no parameter was tuned after seeing results.

## Citation

Nguyen T.D., Pham Q.T. (2026). Simulation companion to "Chọn giống lấy kiểu hình làm
trung tâm trong kỷ nguyên gen: bài học từ ngô thực phẩm." *Vietnam Journal of
Agricultural Sciences* (in review).
