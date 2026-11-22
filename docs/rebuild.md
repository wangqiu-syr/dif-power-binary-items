# Rebuilding the 2011 study in R

This study was first run in 2011, from R through MplusAutomation (Hallquist & Wiley, 2018), with
Mplus 4.2 fitting the MIMIC model and difR's `difLRT()` running the IRT likelihood-ratio test. The original
materials are the author's own code, input and output files, and data. They are kept unchanged in
`original/` of the author's development repository, together with a corrected copy of the 2011
MIMIC input file. The 2011 working folder was named `NCME2011/DIFR_Mpls_MIMIC`.

> **In this repository.** The 2011 materials (`original/` in the author's development repository) and
> the script that reads them, `R/01_single_dataset.R`, are not part of this repository. This page
> describes them and the checks made on them; the check's output is kept in
> [`results/01_single_dataset_output.txt`](../results/01_single_dataset_output.txt) and
> [`results/mimic_lavaan_vs_mplus42.csv`](../results/mimic_lavaan_vs_mplus42.csv).

The rebuild:

- runs entirely in free software (R, lavaan, difR and psychotree);
- is checked against the 2011 output before its results are used;
- fixes the problems found in the original code.

The IRT likelihood-ratio test and the Rasch tree are methods being compared, not cross-checks.
The cross-check is lavaan against the 2011 output (see [Mplus results](#4-mplus-results) below).
The choices made along the way are listed in [decisions.md](decisions.md).

## Step by step

| Step | 2011 original | Issue found | Current version (file / function) | Evidence |
|---|---|---|---|---|
| 1. Data generation | `rgroup()` / `fgroup()` / `respon.1pl()` in `original/Mplus_DIFR_R.r`: 1PL data, 50 + 50 examinees | None in the model. No seed was set (see step 6) | `simulate_dif_data()`, `respon_1pl()` in `R/functions.R` | Identical data to the 2011 code under the same seed, checked over 200 seeds |
| 2. Importing data | `original/DIF_R.txt` reads `dif_responseiteration1.csv`, first with `read.csv(header = TRUE)`, then with `read.table(h = F)` from a hard-coded Windows path (`C:\try1\`) | The files have no header: the first person becomes column names, or each line is read as one column and `data1[, 2:11]` fails | `read_dif_file()` in `R/functions.R`: no header, comma or tab separated | `R/01_single_dataset.R` reads both `dif_responseiteration3.csv` and `DIF_data.dat` ([output](../results/01_single_dataset_output.txt)) |
| 3. Mplus input | `original/diff_ab1.inp`, run by `runModels()` in each replication | Reads the fixed file `dif_responseiteration3.csv`, not `DIF_data.dat`, which the loop rewrote each time: every replication's MIMIC result came from the same data | `fit_mimic()` / `mimic_syntax()` in `R/functions.R`: lavaan, in memory, on each replication's data. `original/diff_ab1_corrected.inp` for re-running with Mplus | `original/diff_ab1.out` shows `Input data file(s): dif_responseiteration3.csv`; the run is in `original/Mplus_Run_Models.log` |
| 4. Mplus results | `original/diff_ab1.out`: Mplus 4.2, 19 Feb 2011, on `dif_responseiteration3.csv` | None; used as the benchmark for lavaan | `R/01_single_dataset.R`, `read_mplus42_results()` in `R/functions.R` | Estimates agree to within 0.001; SEs differ by a factor of about 1.005; same items flagged ([comparison](../results/mimic_lavaan_vs_mplus42.csv)) |
| 5. DIF methods | MIMIC in Mplus (flag `Est./S.E. >= 1.96`) and `difLRT(..., purify = TRUE)` (flag `LRT >= thr`) | The two flag rules use different α levels. MIMIC assumes equal group means on the factor | `run_mimic()`, `run_lrt()`, and the new `run_tree()` in `R/functions.R`; 2011 flag rules kept | Every raw z and LRT statistic is saved, so both can be re-scored at the same α ([replications file](../results/sim_replications_rep1000_MIMIC-LRT-TREE.csv)) |
| 6. Simulation design | `original/Mplus_DIFR_R.r`: 5 replications, no seed | The replications can't be reproduced; 5 replications give flag rates only in steps of 20% | `R/02_simulation.R` (any number of replications, seed set and changeable), `R/03_figures.R` (merges parallel runs) | 1,000 replications, with 95% confidence intervals ([summary](../results/sim_summary_rep1000_MIMIC-LRT-TREE.csv)) |

## 1. Data generation

The data generation is the original 1PL code. The rewritten version, `simulate_dif_data()` and
`respon_1pl()`, reproduces the 2011 `rgroup()` / `fgroup()` draw-for-draw: identical data under
the same seed, checked over 200 seeds. The model and the math are in
[methods.md §2](methods.md#2-the-data-generating-model-rasch-1pl).

The design is the 2011 one:

| | |
|---|---|
| Examinees | 50 reference (`x = 0`) + 50 focal (`x = 1`), θ ~ N(0, 1) in both groups (no impact) |
| Items | 10, Rasch / 1PL |
| Reference b | −1, −1, −1, 0, 0, 0, 0, 1, 1, 1 |
| Focal b | **−1.8**, −1, −1, 0, **−0.8**, 0, 0, 1, 1, **0.2** |
| DIF items | i1, i5, i10: easier for the focal group by 0.8 logits; the other 7 items have no DIF |
| MIMIC flag | `Est./S.E. >= 1.96` for `i ON x` (one-sided, α = .025; the 2011 rule) |
| LRT flag | `difLRT(..., purify = TRUE)`, `LRT >= thr` (χ²(1), α = .05) |
| Tree flag | tree splits on group (α = .05), then anchored Wald test p < .05 (new) |
| Replications | 5 in 2011; 1,000 here |

## 2. Importing data

**`DIF_R.txt` reads its data incorrectly.** It runs difLRT on `dif_responseiteration1.csv`, the
first data set of the 2011 loop. That file is not among the 2011 materials, but the data files of
the same kind (`dif_responseiteration3.csv`, `DIF_data.dat`) have no header.

- `read.csv(..., header = TRUE)` turns the first person into column names.
- `read.table(..., h = F)` has no `sep = ","`, so each line is read as a single column, and
  `data1[, 2:11]` fails.

The second read also uses a hard-coded Windows path, `C:\try1\`. `read_dif_file()` reads both the `.csv` and the
tab-separated `.dat` files.

## 3. Mplus input

**The 2011 MIMIC runs analysed the same data in every replication.** The loop writes each new
data set to `DIF_data.dat`, but `diff_ab1.inp` reads `dif_responseiteration3.csv`. The `.out`
file, written by the run in `Mplus_Run_Models.log`, shows
`Input data file(s): dif_responseiteration3.csv`. So every row of `Pseudo_Z` came from the same
fixed file, and only the LRT results varied across replications. Here, the MIMIC model is fitted
to each replication's own data. `original/diff_ab1_corrected.inp` is a corrected input file for
anyone re-running the original setup.

The flow of the 2011 loop and of this rebuild is drawn in
[methods.md §1](methods.md#2011-pipeline-vs-this-rebuild).

**Lesson.** When an R loop drives an external program through files, check that the program
reads the file the loop just wrote. A one-line check per replication would have caught the 2011
problem: for example, compare the item proportions in each output file with those of the data
just generated. Fitting the model in memory, as lavaan does here, removes the file hand-off
altogether.

## 4. Mplus results

`R/01_single_dataset.R` fits the lavaan model to `dif_responseiteration3.csv` and compares it,
parameter by parameter, with `original/diff_ab1.out`. On this file lavaan reproduces the 2011
MIMIC results:

- **Estimates** (loadings, `i ON x` DIF effects, thresholds) agree to within 0.001.
- **Standard errors** in lavaan are larger by a near-constant factor of about 1.005 ≈ √(N/(N−1)),
  which looks like an N vs. N−1 convention. So z-values are about 0.5% smaller; for example, for
  i5, z = 2.413 in the 2011 output and 2.401 in lavaan.
- **DIF decisions** are the same: i5 and i10 are flagged in both. i1 is not flagged
  (z = 1.51 in the 2011 output, 1.50 in lavaan).
- **Model fit.** The 2011 output reports the older WLSMV χ² (mean- and variance-adjusted,
  estimated df), which lavaan reproduces closely. lavaan's WLSMV default is the newer
  scaled-and-shifted test.

| | 2011 output | lavaan, older test | lavaan WLSMV default |
|---|---|---|---|
| χ² | 29.252 | 29.163 | 43.809 |
| df | 22 | 22.155 | 35 |
| p | .138 | .145 | .146 |
| CFI | .939 | .938 | .946 |
| TLI | .928 | .940 | .931 |
| RMSEA | .057 | .057 | .050 |
| WRMR | .743 | .739 | .739 |

## 5. DIF methods

**MIMIC** is a one-factor model for the 10 binary items, with each item regressed on group
(`i ON x`), estimated with WLSMV. It is fitted in lavaan (Rosseel, 2012), with the same model and
estimator as the 2011 Mplus run. The **IRT likelihood-ratio test** is `difR::difLRT()` (Magis et
al., 2010) with purification, called as in the 2011 code ([difLRT.md](difLRT.md)). The **Rasch
tree** is new: `psychotree::raschtree()` (Strobl, Kopf, & Zeileis, 2015) tests whether the Rasch
item parameters differ between groups, and if they do, anchored Wald tests from
`psychotools::anchortest()` (Kopf, Zeileis, & Strobl, 2015) pick out the DIF items. Tree methods
connect to the author's earlier work on classification trees (Vaughn & Wang, 2008;
[methods.md §7](methods.md#7-how-the-rasch-tree-relates-to-the-2008-classification-trees)).

**The two flag rules use different α levels.** The MIMIC rule `z >= 1.96` is one-sided
(α = .025) and only catches DIF that favours the focal group. The LRT rule is a two-sided χ²(1)
test at α = .05. The 2011 rules are kept for comparability. All raw statistics are saved, so both
methods can be re-scored at the same α.

**A modelling assumption to keep in mind.** The MIMIC model frees all ten direct effects and has
no `f1 ON x`, which identifies it by assuming equal group means on the factor. That holds in this
design. If impact is added, add `f1 ~ x` and fix at least one item's direct effect to zero as an
anchor. The LRT allows a group mean difference either way.

## 6. Simulation design

**No seed was set.** The 2011 replications can't be reproduced. A seed is now set and can be
changed. **Five replications** give flag rates only in steps of 20%. Any number can now be run.

The Monte Carlo approach follows Vaughn and Wang (2008): simulate in R, verify the code against
known results before simulating, and use enough replications for stable error rates.

The reported run: 1,000 replications, run as 8 parallel jobs of 125 replications
(`R/02_simulation.R 125 <seed>`, seeds 101–108) and merged by `R/03_figures.R`. All three methods
analyse the same 1,000 data sets. R 4.5.1, lavaan 0.7.2, difR 6.1.0, psychotree 0.16.3,
psychotools 0.7.7. All fits ran; 9 MIMIC fits gave lavaan warnings and are kept.

Flag rate per item, from [`sim_item_rates_rep1000_MIMIC-LRT-TREE.csv`](../results/sim_item_rates_rep1000_MIMIC-LRT-TREE.csv)
(drawn in [Figure 2](../results/figures/fig2_item_flag_rates.png)):

| item | i1 | i2 | i3 | i4 | i5 | i6 | i7 | i8 | i9 | i10 |
|---|---|---|---|---|---|---|---|---|---|---|
| DIF | yes | | | | yes | | | | | yes |
| MIMIC | .292 | .024 | .025 | .032 | .385 | .026 | .033 | .021 | .022 | .366 |
| IRT-LRT | .283 | .115 | .132 | .122 | .346 | .119 | .125 | .122 | .120 | .339 |
| Rasch tree | .106 | .029 | .053 | .042 | .132 | .047 | .043 | .044 | .032 | .126 |

The Rasch tree split on group in 28.7% of replications. Among those, its power was .42 and its
Type I error .14.

## References

Hallquist, M. N., & Wiley, J. F. (2018). MplusAutomation: An R package for facilitating large-scale
latent variable analyses in Mplus. *Structural Equation Modeling, 25*(4), 621–638.

Kopf, J., Zeileis, A., & Strobl, C. (2015). A framework for anchor methods and an iterative forward
approach for DIF detection. *Applied Psychological Measurement, 39*(2), 83–103.

Magis, D., Béland, S., Tuerlinckx, F., & De Boeck, P. (2010). A general framework and an R package
for the detection of dichotomous differential item functioning. *Behavior Research Methods, 42*(3), 847–862.

Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
*Journal of Statistical Software, 48*(2), 1–36.

Strobl, C., Kopf, J., & Zeileis, A. (2015). Rasch trees: A new method for detecting differential
item functioning in the Rasch model. *Psychometrika, 80*(2), 289–316.

Vaughn, B. K., & Wang, Q. (2008). Classification based on tree-structured allocation rules.
*The Journal of Experimental Education.*
