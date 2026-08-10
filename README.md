# Simulation study: major-gene-fixation variance depletion in edible maize breeding

Companion simulation study for the review "Chọn giống lấy kiểu hình làm trung tâm trong
kỷ nguyên gen: bài học từ ngô thực phẩm" (*Vietnam Journal of Agricultural Sciences*).
It turns the review's central quantitative-genetics claim — that once a large-effect
quality locus (sh2-like) is fixed in a commercial market class, the remaining additive
genetic variance sits entirely in the polygenic background, and different post-fixation
selection strategies exploit that residual variance with different efficiency — from a
narrative argument into a reproducible, statistically replicated demonstration.

A full Nature-style manuscript write-up (Methods, Results, Discussion) is in
[`MANUSCRIPT.md`](MANUSCRIPT.md). This README summarizes the same study.

Built with [AlphaSimR](https://cran.r-project.org/package=AlphaSimR) (Gaynor et al.
2021, *The Plant Genome*/*G3*), [AGHmatrix](https://cran.r-project.org/package=AGHmatrix),
and [sommer](https://cran.r-project.org/package=sommer).

## What it simulates

**Phase 1 — historical fixation.** A founder population (N = 800, 10 chromosomes,
coalescent LD via `runMacs`) segregating at one large-effect quality locus plus 300
small-effect background QTL. Eighteen generations of truncation selection (top 20%) on
the composite phenotype fix the major locus, at which point its own additive variance
collapses to zero — matching the classical prediction V_A = 2pq·a² closely at every
generation. Background additive variance, by contrast, stays firmly above zero — the
literal, numeric version of the review's "cạn biến dị do cố định gen lớn" argument.

**Phase 2 — three post-fixation strategies**, all starting from Phase 1's final
(major-locus-fixed) population and run for 12 further generations:

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

## Why replication, not a single run

AlphaSimR's founder-haplotype generator (`runMacs`) calls an external coalescent sampler
whose internal RNG state is not fully pinned by R's own `set.seed()`. A single run is
therefore not exactly reproducible, and — more importantly — is not on its own a
statistically defensible basis for comparing three strategies. The full Phase 1 → Phase
2 pipeline was independently repeated **20 times**, and every number below is a mean ±
SD across those 20 independent replicates, with paired t-tests (each replicate
contributes one paired observation per arm, since all three arms in a replicate share
the same Phase 1 starting population).

## Real results (read directly from `outputs/*.csv`, n = 20 replicates)

**Phase 1 — fixation and variance collapse:**

| Quantity | Value |
|---|---:|
| Generation to fixation (p ≥ 0.99), mean ± SD | 3.4 ± 1.4 |
| Additive variance at major locus, initial mean | 1.001 |
| Additive variance at major locus, at fixation | 0.0000 |
| Background additive variance, at major-locus fixation | 0.392 |
| Background additive variance, after 18 total generations | 0.211 |

**Phase 2 — three strategies:**

| Strategy | Cumulative genetic gain, mean ± SD | Relative to A |
|---|---:|---:|
| A — Phenotypic mass selection | 3.801 ± 0.308 | 1.00× |
| B — Multi-trait index selection | **4.063 ± 0.345** | **1.07×** |
| C — Genomic selection (GBLUP) | 3.046 ± 0.276 | 0.80× |

Multi-trait index selection (B) beat phenotypic mass selection (A) by **6.9%**
(paired t-test, t(19) = 10.10, p < 0.001) and beat genomic selection (C) by **33.4%**
(paired t-test, t(19) = 28.29, p < 0.001). Mean genomic-selection prediction accuracy
across cycles was **r = 0.319 ± 0.020** — modest, because the training population here
is freshly re-sampled each cycle (no historical accumulation) and the marker panel is
moderate relative to the 300-locus background architecture. This is a realistic
constraint for a resource-limited program, not a design flaw: it is exactly why genomic
selection's real published advantage (e.g. Beyene et al. 2015, cited in the main review)
comes mainly from *faster cycling*, not higher per-cycle accuracy — a nuance this
simulation reproduces on its own rather than assuming.

**Net finding**: once the major locus is fixed, *how you structure phenotyping* —
measuring the true quality sub-components separately and combining them with a proper
index — recovered more of the residual polygenic variance than either naively
phenotyping one composite score, or defaulting to marker-based prediction with an
under-resourced training population. This is the review's own Section 7 argument, now
with a statistically replicated number attached.

## Reproducing

```r
install.packages(c("AlphaSimR", "AGHmatrix", "sommer", "ggplot2", "dplyr", "patchwork"))
setwd("scripts")
source("01_replicated_pipeline.R")  # runs 20 independent replicates; writes outputs/*_replicated_*.csv
source("02_nature_figures.R")       # writes figures/Figure1_*, figures/Figure2_*
source("03_tables.R")               # writes outputs/tables/Table1-3*.csv, prints MANUSCRIPT.md values
```

Each replicate reseeds explicitly (`set.seed(20260810 + rep_id)`), but because
`runMacs()`'s own coalescent sampler is not fully controlled by R's seed (see above),
individual replicate values are **not** expected to reproduce bit-for-bit on a re-run —
the aggregate statistics across 20 replicates are the reported, stable result, and are
consistent with the numbers above.

## Caveats

- Sub-trait weights in Arm B use the *true* simulated genetic variances (a known-parameter
  index), isolating the value of decomposition + optimal weighting from the separate,
  real-world problem of estimating those weights from limited data.
- Arm C retrains from a fresh random half of the population each cycle; a program that
  accumulates training data across cycles, or that increases marker density, would likely
  close some or all of the gap with B — this simulation shows one realistic operating
  point, not a universal ranking of methods.
- Population size (N = 800), generation counts, and heritabilities were fixed a priori
  before any replicate was run; no parameter was tuned after seeing results.
- Reported allele frequencies at the major locus are direction-agnostic
  (`max(p, 1-p)`, i.e. "distance to fixation"), because `AlphaSimR::addTraitA()` assigns
  the effect sign to the reference allele arbitrarily per replicate — raw frequencies are
  not directly comparable/averageable across replicates without this transform.

## Citation

Nguyen T.D., Pham Q.T. (2026). Simulation study: major-gene-fixation variance depletion
in edible maize breeding — companion to "Chọn giống lấy kiểu hình làm trung tâm trong
kỷ nguyên gen: bài học từ ngô thực phẩm." *Vietnam Journal of Agricultural Sciences* (in
review).
