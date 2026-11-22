# Decisions

The choices made in rebuilding the study, as described in [README](../README.md),
[methods.md](methods.md), [difLRT.md](difLRT.md) and [rebuild.md](rebuild.md). Impact: how much
the choice affects the results or their interpretation (H = high, M = medium, L = low).

| ID | Topic | Choice | Impact |
|---|---|---|---|
| D1 | MIMIC software | lavaan instead of Mplus | H |
| D2 | Flag rules | Keep the 2011 rules: MIMIC `z >= 1.96` (one-sided, α = .025), LRT `LRT >= 3.84` (two-sided, α = .05) | H |
| D3 | MIMIC identification | No `f1 ON x`; all ten direct effects free; equal group means assumed | M |
| D4 | IRT-LRT package | difR's `difLRT()`, as in 2011 | M |
| D5 | difLRT purification | `purify = TRUE`, as in 2011 | M |
| D6 | Third method | Add the Rasch tree (psychotree), with `anchortest()` (forward class, MTT selection) for item-level tests | M |
| D7 | Data generation | Keep the original 1PL code, reproduced draw-for-draw | M |
| D8 | Seeds and random numbers | Set a seed; each data set depends only on the seed and the replication number; per-replication seed for the tree's tie-breaking | M |
| D9 | Replications | 1,000 instead of 5, run in parallel and merged | M |
| D10 | Failed fits and warnings | Failed fits (NA) are left out of the denominator; fits with lavaan warnings are kept | L |
| D11 | The 2011 materials | Kept unchanged in `original/` of the development repository (not in this repository); the corrected input file is a separate copy | L |
| D12 | Power ceilings | Computed exactly from the 1,024 response patterns, not by simulation | M |
| D13 | difLRT for the diagnosis | Rebuilt in lme4 from difR 6.1.0's source, with one optimizer (bobyqa) | M |
| D14 | Power pilot | Same data sets as the 1,000-replication run (seeds 101–102); 100 replications, 200 at n = 150 | M |
| D15 | A correctly computed LRT | Conditional ML (psychotools) instead of glmer deviance differences | H |
| D16 | Mantel–Haenszel in the pilot | With and without the continuity correction | L |
| D17 | Power study size | 40 cells × 300 replications (planned 500); difLRT in 2 cells only | M |
| D18 | Robust anchor | The k items closest to the median CML group difference; k = 5 and k = 3 both run | H |
| D19 | Flag rules in the power study | Two-sided .05 primary; one-sided .05, Benjamini–Hochberg and size-adjusted power from saved statistics | M |
| D20 | Storing the power study's statistics | Compact `replications.rds` (statistics to 4 decimals, p to 4 significant digits) instead of a 24 MB CSV | L |
| D21 | Mantel–Haenszel variance in double precision | `mh_one()` computed nR·nF·m1·m0 in integers, which overflow above about 430 examinees per score level; now double. No effect on the power study (largest strata well below that; 0 missing MH results in 120,000) | L |
| D22 | CML fits with fewer than four items per missing-data pattern | Numeric-Hessian refit only when psychotools 0.7-3's analytic Hessian fails; no effect on the power study | L |

## D1. lavaan instead of Mplus

The 2011 study fitted the MIMIC model in Mplus 4.2. The rebuild uses lavaan, with the same
specification and the same WLSMV estimator, verified against the 2011 output to within 0.001
([rebuild.md §4](rebuild.md#4-mplus-results)). lavaan runs inside the R loop, with no data files
to write, external program to call, or output to parse, which removes the step where the 2011
bug happened. It is free, so anyone can rerun the study without an Mplus licence
([methods.md §6](methods.md#6-why-difr-lavaan-and-psychotree)).

## D2. Keep the 2011 flag rules

The MIMIC rule `z >= 1.96` is one-sided (α = .025) and only catches DIF that favours the focal
group; the LRT rule is a two-sided χ²(1) test at α = .05. The rules are kept for comparability
with 2011. At face value, MIMIC is held to a stricter false-positive rate. All raw z and LRT
statistics are saved, so both methods can be re-scored at the same α without rerunning.

## D3. MIMIC identification

The MIMIC model frees all ten direct effects and has no `f1 ON x`, which identifies it by
assuming equal group means on the factor. That holds in this design, where both groups are
N(0, 1). If impact is added, add `f1 ~ x` and fix at least one item's direct effect to zero as
an anchor. The LRT allows a group mean difference either way.

## D4. difR for the IRT likelihood-ratio test

difR is the peer-reviewed R package for dichotomous DIF (Magis et al., 2010). The IRT
likelihood-ratio test is the method MIMIC is usually compared with, its Rasch model matches the
data-generating model, and the 2011 code used `difLRT()`, so keeping it keeps the replication
faithful ([methods.md §6](methods.md#6-why-difr-lavaan-and-psychotree)).

## D5. difLRT purification

`purify = TRUE`, as in the 2011 code. Flagged items are dropped from the anchor set and the tests
are repeated until the flagged set stops changing, or until `nrIter` iterations (default 10).
This reduces the distortion that DIF in anchor items causes
([difLRT.md](difLRT.md#purification-purify--true)). Purification makes difLRT much slower than
lavaan.

## D6. Add the Rasch tree

The rebuild adds `psychotree::raschtree()` as a third method. If the tree splits on group,
`psychotools::anchortest()` with the iterative forward class and MTT selection picks the anchor
and Wald-tests each item. The first anchor item fixes the scale, so it gets no test and is
counted as not flagged; if the tree doesn't split, no item is flagged
([methods.md §5](methods.md#5-rasch-tree-psychotree-added-in-this-rebuild)). psychotree is the
reference implementation by the method's authors.

## D7. Keep the original data generation

`simulate_dif_data()` and `respon_1pl()` reproduce the 2011 `rgroup()` / `fgroup()` draw-for-draw:
identical data under the same seed, checked over 200 seeds
([rebuild.md §1](rebuild.md#1-data-generation)).

## D8. Seeds and random numbers

The 2011 code set no seed, so its replications can't be reproduced. `R/02_simulation.R` sets a
seed (default 20110219), which can be changed. The random-number state is saved after each data
set is generated and restored after the methods have run, so each data set depends only on the
seed and the replication number, not on which methods are run. Runs of different methods with
the same seed analyse the same data sets. `anchortest()` breaks ties at random, so the tree gets
its own seed in each replication (`seed + 100000 + r`).

## D9. 1,000 replications, in parallel

Five replications give flag rates only in steps of 20%. The reported run has 1,000, run as 8
parallel jobs of 125 (seeds 101–108) and merged by `R/03_figures.R`
([rebuild.md §6](rebuild.md#6-simulation-design)).

## D10. Failed fits and warnings

A failed or non-converged fit returns NA instead of stopping the simulation, and NAs are left out
of the denominator of the flag rates. MIMIC fits that give lavaan warnings are kept; in the
reported run, 9 of 1,000 did.

## D11. The 2011 materials stay unchanged

The original code, input and output files, and data are kept unchanged in `original/` of the
author's development repository; they are not part of this repository ([rebuild.md](rebuild.md)).
The corrected MIMIC input, which reads the data file the loop writes, is a separate file,
`original/diff_ab1_corrected.inp`, for anyone re-running the original setup.

## D12. Power ceilings computed exactly

With 10 binary items there are 1,024 response patterns, so each group's pattern probabilities can
be computed by numerical integration over θ, and a model can be fitted to them directly. Twice
the gain in expected log-likelihood is the noncentrality of the LRT, so the expected power of
each kind of test, and the sample size it needs, follow without simulation and without Monte Carlo
error ([power.md §2](power.md#2-how-much-power-can-50--50-examinees-give); results in
`results/power/exact_*.csv`, code released with the follow-up study).
The pilot checks these numbers.

## D13. difLRT rebuilt in lme4 for the diagnosis

difR could not be installed on the Linux machine used for the diagnosis (no access to CRAN), so `LRT()` and the
purification loop of `difLRT()` were rebuilt from difR 6.1.0's source (`R/power/lrt_variants.R`).
The rebuild reproduces difR's first-round statistics to four decimals with lme4's default
optimizer. The pilot uses bobyqa alone, which is about 4 times faster; on the 100 shared data sets
its flags agree with the laptop's difR run in 99.8% of item tests (correlation of the statistics
.9997).

## D14. The power pilot

The pilot reuses the first 50 data sets of seeds 101 and 102 of the 1,000-replication run, so its
results can be read against the laptop's (tree splits and MIMIC z are identical). It adds a design
with no DIF (seed 301, for Type I error and size-adjusted power) and the 2011 DIF at 150 per group
(seed 401). 100 replications per design (200 at n = 150) give power to within about ±.05 and Type
I error to within about ±.02; that is enough to locate the losses, not to rank close methods.

## D15. A correctly computed LRT

difR's `LRT()` takes the difference of lme4's `deviance()` between the two models. For a GLMM,
lme4's `deviance()` is the conditional deviance at the predicted random effects, not −2
log-likelihood, so the difference is not a likelihood ratio. With no DIF its mean is 1.17–1.26
instead of 1, and it rejects 6–8% of the time instead of 5%
([power.md §5](power.md#5-pilot-where-the-power-goes)). The improved LRT therefore uses conditional
maximum likelihood (`psychotools::raschmodel`): the studied item is split into a reference and a
focal "virtual item", and the LR statistic compares the two fits. CML needs no assumption about
the θ distributions, so impact needs no model, and each fit takes milliseconds.

## D16. Mantel–Haenszel with and without the continuity correction

difR's `difMH()` applies the continuity correction by default. At 50 + 50 it removes 0.5 from a
difference whose standard deviation is about 2, a large share, so the pilot reports both.

## D17. Power study size

The plan in [power.md §7](power.md#7-proposed-full-study) had 500 replications per cell.
Timing in the pilot (about 4 s per replication at 50 per group and 12 s at 250, mostly MIMIC with
scale purification) put 500 at about 13 hours on 4 cores, so the study runs 300 (about 7–8
hours): power to within about ±.03 and Type I error to within about ±.01 per cell. difLRT (glmer,
about 10 s more per replication) runs only in the 2011 design and its no-DIF twin, where it is the
reference; its behaviour elsewhere follows from power.md §3.2. MIMIC fits use `test = "none"`,
which leaves the direct effects' z unchanged and saves about a third of the time.

## D18. The robust anchor

Each group's Rasch difficulties are estimated by CML and centred; d_i = b_iR − b_iF. The anchor for
item j is the k items (other than j) closest to the median of the d_i, and item j is tested with
the CML LRT against it ([power_study.md §3](power_study.md#3-methods-each-run-on-every-data-set)).
The median is not moved by DIF in fewer than half the items, which is the idea behind robust
scaling (Halpin, 2024) and the L1 condition of Chen et al. (2023); this is a simple version of it,
not those authors' estimators. k trades power (a longer anchor) against the risk of including a
DIF item. Both k = 5 and k = 3 are run and reported; neither was tuned on the results.

## D19. Flag rules in the power study

Two-sided α = .05 is the primary rule, the same for every method. Because every signed statistic
and p-value is saved, the one-sided rule (hypothesized direction; one-way and no-DIF cells only),
Benjamini–Hochberg control across the 10 items, and size-adjusted power (critical value from the
no-DIF cell with the same n and impact) are applied afterwards, without rerunning.

## D20. Storing the power study's statistics

The full table of statistics (1.1 million rows) is 24 MB as a gzipped CSV with the design columns
repeated on every row. `R/09_power_study_summary.R` instead writes `replications.rds` (cell, rep,
method, item, statistic rounded to 4 decimals, p-value to 4 significant digits; 4.7 MB with xz
compression) and a per-replication `replications_info.csv`. The design columns come from
`study_cells()` by cell. The rates in `rates.csv` are computed from the unrounded statistics and
are identical either way (checked). The per-replication `.rds` files of the run stay untracked;
rerunning `R/08_power_study.R` with the same seeds regenerates them.

## D21. Mantel–Haenszel variance in double precision

`mh_one()` in `R/power/pilot_methods.R` summed the hypergeometric variance term
`nR * nF * m1 * m0` from integer counts. R integers overflow above 2^31 − 1, which this product
reaches once a score level holds more than about 430 examinees, as in large survey samples. The term is now
computed in double precision. Results are identical where there is no overflow (checked on a
50 + 50 data set to six decimals), and the power study never reached it: its largest score levels
hold well under 430 examinees, and none of its 120,000 Mantel–Haenszel results is missing.

## D22. CML fits with fewer than four items per missing-data pattern

psychotools 0.7-3 (the version used throughout, R 4.3.3) fails in `raschmodel()` when a
missing-data pattern has fewer than four items: its analytic Hessian calls
`elementary_symmetric_functions()` with one parameter left and stops with "n > 1L is not TRUE".
This happens in a CML likelihood-ratio test once purification leaves an anchor of two items (two
anchor items plus the split item per pattern). A later application to survey data reached it: standard purification
flagged six of eight items, and every test of the next pass returned NA. `rasch_cml()` in
`R/power/study_methods.R` now refits with the numeric Hessian only when that error occurs. The
estimates and log-likelihood are the same (the Hessian is not used for them); the standard errors
agree with the analytic ones to five decimals where both run (checked on a 4-item data set). Fits
that worked before take the unchanged path, so earlier results cannot change; the power study
never reached the error (no CML purified test returned NA, `results/power_study/fit_problems.csv`).
