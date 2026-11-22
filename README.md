# A wider window, the same wall: Two decades of item-level DIF tests for binary items

MIMIC, IRT likelihood-ratio and Rasch-tree tests for binary items: a report, with the code and
results behind every number in it.

[Qiu Wang](https://orcid.org/0000-0003-1635-8050) (Syracuse University) · Version 1.0.1 ·
November 2026 · DOI (version 1.0.0): [10.5281/zenodo.24447165](https://doi.org/10.5281/zenodo.24447165) · Preprint: to be added

- **Report:** [`report/report.pdf`](report/report.pdf) (source [`report/report.md`](report/report.md))
- **Supplement A3**, the retrospective ADEMP-PreReg form (Siepe et al., 2024) for Simulation 2:
  [`report/supplement_A3.pdf`](report/supplement_A3.pdf) (source [`report/supplement_A3.md`](report/supplement_A3.md))

## Key findings

1. **The baseline: at 50 examinees per group, every method has low power.** With 10 Rasch items and
   0.8 logits of DIF on 3 of them, MIMIC flagged a DIF item 34.8% of the time, the IRT
   likelihood-ratio test (difLRT) 32.3% and the Rasch tree 12.1% (1,000 replications). Exact
   calculations show that no valid two-sided test of a single item can exceed about .39 there: the ceiling is the design.
2. **Three methods become nine.** The three methods differ in four components at once (statistic,
   anchor, design, decision rule). A pre-specified simulation (40 conditions × 300 replications)
   used nine methods in three families (LRT, MIMIC, Rasch tree), each variant changing one or two
   components, plus a θ-known benchmark, so that each factor could be isolated. Table 11 rates every
   method on the twelve factors (Yes / Semi / No).
3. **What raises power.** Sample size above all: for 0.8 logits, .80 power per item needs about
   160–195 examinees per group with a purified anchor (about 150–160 with a known clean anchor);
   for 0.5 logits, even a known clean anchor needs about 380. At a fixed sample size, a one-sided
   test in a pre-registered direction (+.07 to +.12) and a joint test of the suspected items (.70
   at 50 per group) raise power; a median-based robust anchor did not.
4. **What controls Type I error.** A correctly computed statistic, purification that re-tests every
   item and, for MIMIC, modeling impact. Purified CML LRT, purified Mantel–Haenszel and MIMIC with
   scale purification are nearly interchangeable and hold about .05 without DIF or with mixed DIF
   (.05–.09 with one-way DIF); difLRT's extra power was false
   positives, and the baseline MIMIC flagged 17–60% of clean items when the groups differed by
   0.5 SD.

## Overview

The code in this repository reproduces every table and figure of the report from two Monte Carlo
simulations, both on simulated data. This README follows the template of Vilhuber et al. (2022).

- **Simulation 1** (report §3): the baseline design, 50 + 50 examinees, 1,000 replications,
  three methods (MIMIC in lavaan, `difR::difLRT`, the Rasch tree in psychotree).
- **Simulation 2** (report §5): the pre-specified power study, 40 conditions × 300 replications
  (12,000 data sets), nine methods in three families plus a θ-known benchmark, a three-item robust
  anchor and a bundle test, hypotheses H1–H6 fixed before the results.
- In between (report §4): a pilot of the remedies (`R/04_power_pilot.R`) and a calibration check
  of difLRT's statistic (`R/07_lrt_calibration.R`).

`Rscript R/run_all.R` rebuilds the report's tables and figures from the stored results and runs
the tests in under a minute. `Rscript R/run_all.R FULL` reruns both simulations (about 12 hours on
8 cores).

## Data availability and provenance

All data are simulated by the code in `R/`; no external data are used, and there are no access
restrictions. The stored results are in `results/`. Every simulated data set can be regenerated
from its seed (see *Controlled randomness*).

### What is not in this repository

- **The code for the exact power calculations** (`R/power/exact_power.R`, `R/05_power_exact.R`).
  Its results, which the report uses (Table 4, §4.1, Table 5, Figure 3), are in
  `results/power/exact_*.csv`. The code is released with the follow-up study, registered on OSF
  (https://osf.io/35nvh/).
- **The original materials of the baseline study** (Mplus and difR code, Mplus input and output,
  data files) and the script that compared lavaan with the original Mplus output
  (`R/01_single_dataset.R`).
  [`docs/rebuild.md`](docs/rebuild.md) describes them and the checks made; the check's output is
  kept in `results/01_single_dataset_output.txt` and `results/mimic_lavaan_vs_mplus42.csv`.

That is why the script numbers skip 01 and 05.

## Computational requirements

### Software

- **R 4.3.3** with lavaan 0.6-17, psychotools 0.7-3, psychotree 0.16-0, partykit 1.2-20,
  mvtnorm 1.2-4 and lme4 1.1-35.1, on Linux (Ubuntu 24.04): the versions of the pilot, the
  calibration check and Simulation 2
  ([`results/power/sessionInfo.txt`](results/power/sessionInfo.txt),
  [`results/power_study/sessionInfo.txt`](results/power_study/sessionInfo.txt)). These runs compute
  difLRT's statistic with an lme4 rebuild of difR 6.1.0 (`R/power/lrt_variants.R`), so they do not
  need difR.
- Simulation 1 was run on a Windows laptop with R 4.5.1 (lavaan 0.7.2, difR 6.1.0, psychotree
  0.16.3). The tests check that R 4.3.3 on Linux reproduces its MIMIC and Rasch-tree statistics
  (and its difLRT statistics when difR is installed).
- `R/report_tables.R`, `R/report_figures.R` and the tests' report checks use base R only.
- Building the PDF needs pandoc (3.1 or later), Chromium or Google Chrome, and an OpenType math
  font for the formulas (TeX Gyre Termes Math, Ubuntu package `fonts-texgyre-math`; STIX Two Math
  or Cambria Math also work).

Install the packages once, from R:

```r
install.packages(c("lavaan", "psychotools", "psychotree", "partykit", "mvtnorm", "lme4", "difR"))
```

### Controlled randomness

Every script sets its own seeds; nothing depends on the order in which parallel jobs finish.

| Part | Script | Seeds |
|---|---|---|
| Simulation 1 | `R/02_simulation.R` | 8 jobs of 125 replications, seeds 101–108 (`set.seed(seed)` once, then the r-th data set of the stream); the Rasch tree's tie-breaking: `seed + 100000 + r` |
| Pilot | `R/04_power_pilot.R` | Baseline design: seeds 101 and 102 (the same data sets as Simulation 1); no-DIF design: 301; 150 per group: 401; per job `seed × 1000 + r` |
| LRT calibration | `R/07_lrt_calibration.R` | 501 |
| Simulation 2 | `R/08_power_study.R` | replication r of condition c: `20261119 + 1000c + r`; anchortest's tie-breaking: that seed + 7 |
| Paired MCSEs | `R/11_paired_mcse.R` | none (reads the stored statistics) |
| Robust-anchor check | `R/10_robust_anchor_check.R` | Simulation 2's seeds for condition 25; `set.seed(2)` for the normal approximation |
| Tests | `tests/` | the seeds above, to regenerate stored replications |

### Memory, run time and storage

| Step | Command | Run time |
|---|---|---|
| QUICK: tables, figures, tests | `Rscript R/run_all.R` | under a minute |
| Simulation 1 | `Rscript R/run_all.R FULL 8 sim1` | about 4 hours (8 jobs in parallel; about a minute per replication, mostly difLRT) |
| Pilot | `Rscript R/run_all.R FULL 4 pilot` | not recorded; about 1–2 hours on 4 cores (estimate) |
| LRT calibration | `Rscript R/run_all.R FULL 4 calibration` | about 10 minutes on 4 cores |
| Simulation 2 | `Rscript R/run_all.R FULL 4 study` | about 7 hours on 4 cores (27.5 CPU hours; 8.2 s per replication) |

Memory use was not measured; the runs used a laptop and a 4-core Linux machine. The stored results
take 9 MB. A FULL run also writes one file per replication to `results/power/runs/` and
`results/power/study_runs/` (not tracked; tens of MB).

## Instructions

From the repository root:

```bash
Rscript R/run_all.R
```

QUICK mode rebuilds Tables 3, 4, 8, 9 and 11 (`results/report_tables.md`) and Figures 1–5
(`report/figures/`) from the stored results, and runs `tests/run_tests.R`. The tests regenerate
replications 1–2 of Simulation 1 (seed 101) and replication 1 of Simulation 2's condition 1 and
compare them with the stored statistics, and check that every generated table row appears in the
report.

To rerun everything, overwriting the results:

```bash
Rscript R/run_all.R FULL 8
```

The steps can also be run one at a time (`Rscript R/run_all.R FULL 8 sim1,study`). An interrupted
run can be restarted: finished replications are skipped. On Windows the eight Simulation 1 jobs
run one after another.

Then build the PDFs (see *Build the PDF*).

## List of tables and figures

| Exhibit | Content | Script | Output / source file |
|---|---|---|---|
| Table 1 | The three baseline methods, by component | — (text) | — |
| Table 2 | ADEMP summary of Simulation 1 | — (text) | — |
| Table 3 | Power and Type I error, Simulation 1 | `R/02_simulation.R`, `R/03_figures.R` → `R/report_tables.R` | `results/sim_summary_rep1000_MIMIC-LRT-TREE.csv` |
| Table 4 | Exact power at 50 + 50, n for .80 | exact calculations (code with the follow-up study) → `R/report_tables.R` | `results/power/exact_items.csv`, `exact_n_for_power.csv` |
| Table 5 | Factors that raise power and factors that control Type I error | — (text, from the exact results) | `results/power/exact_*.csv`, `lrt_calibration.csv` |
| Table 6 | ADEMP summary of Simulation 2 | — (text) | Supplement A3 |
| Table 7 | The nine methods by family (LRT, MIMIC, Rasch tree): what each changes, the factor it tests | — (text) | — |
| Table 8 | Type I error, Simulation 2 | `R/08_power_study.R`, `R/09_power_study_summary.R` → `R/report_tables.R` | `results/power_study/rates_mcse.csv` |
| Table 9 | Power, Simulation 2 | as Table 8 | `results/power_study/rates_mcse.csv` |
| Table 10 | Effect of each factor | `R/09_power_study_summary.R`, `R/10_robust_anchor_check.R` | `results/power_study/rates.csv`, `robust_anchor_check.txt` |
| Table 11 | Support of each method for the twelve factors (Yes / Semi / No) | `R/report_tables.R` | `results/power_study/factor_matrix.csv` (metrics and ratings) |
| Table 12 | Recommendations | — (text) | — |
| Figure 1 | Power and Type I error, Simulation 1 | `R/report_figures.R` | `report/figures/figure1.png` from `results/sim_summary_rep1000_MIMIC-LRT-TREE.csv` |
| Figure 2 | Flag rate per item, Simulation 1 | `R/report_figures.R` | `report/figures/figure2.png` from `results/sim_item_rates_rep1000_MIMIC-LRT-TREE.csv` |
| Figure 3 | Pilot: simulated against predicted | `R/report_figures.R` | `report/figures/figure3.png` from `results/power/pilot_rates.csv`, `exact_items.csv`, `exact_contamination_by_n.csv` |
| Figure 4 | Type I error, Simulation 2 | `R/report_figures.R` | `report/figures/figure4.png` from `results/power_study/rates.csv` |
| Figure 5 | Power, Simulation 2 | `R/report_figures.R` | `report/figures/figure5.png` from `results/power_study/rates.csv` |
| Appendix A1–A2 | All conditions with MCSEs; paired MCSEs of differences between methods | `R/report_tables.R`, `R/11_paired_mcse.R` | `results/power_study/rates_mcse.csv`, `paired_diff_mcse.csv` |

The simulation scripts also draw working versions of the figures in `results/figures/`
(`fig1`–`fig5`); the report uses the redrawn ones in `report/figures/`.

## Repository layout

| Path | Contents |
|---|---|
| `report/` | `report.md`, `supplement_A3.md`, their PDFs, `figures/`, and `pdf/` (the PDF build) |
| `R/functions.R` | Data generation; MIMIC (lavaan), difLRT and Rasch-tree wrappers |
| `R/02_simulation.R`, `R/03_figures.R` | Simulation 1; merging the parallel runs, pooled rates with 95% CIs, working figures |
| `R/04_power_pilot.R`, `R/06_power_summary.R`, `R/07_lrt_calibration.R` | The pilot and its summaries; the calibration of difLRT's statistic |
| `R/08_power_study.R`, `R/09_power_study_summary.R`, `R/10_robust_anchor_check.R`, `R/11_paired_mcse.R` | Simulation 2: the run, the summaries (every flag rule), the robust-anchor diagnostic, paired MCSEs |
| `R/power/` | `lrt_variants.R` (difLRT rebuilt in lme4), `pilot_methods.R` (benchmarks, Mantel–Haenszel, MIMIC variants, the tree with and without its gate), `study_methods.R` (CML LRTs, robust anchor, bundle test, M-SP) |
| `R/report_tables.R`, `R/report_figures.R` | The report's tables (with MCSEs) and figures, from `results/` only |
| `R/run_all.R` | QUICK / FULL entry point |
| `results/` | Simulation 1 (`sim_*_rep1000_MIMIC-LRT-TREE.csv`), `power/` (exact results, pilot, calibration), `power_study/` (Simulation 2: rates, MCSEs, every statistic in `replications.rds`, the pilot), `report_tables.md`, working figures |
| `docs/` | `methods.md` (models and code), `difLRT.md`, `power.md` (why power is low; the plan with H1–H5), `power_study.md` (Simulation 2 in full), `rebuild.md` (the rebuild of the baseline study), `decisions.md` (D1–D22) |
| `tests/` | `run_tests.R` and the tests |

## Build the PDF

```bash
bash report/pdf/build.sh
```

This converts `report/report.md` and `report/supplement_A3.md` to HTML with pandoc
(`report/pdf/filter.lua`, `report/pdf/style.css`; formulas as MathML) and prints them to PDF with headless Chromium
(`report/pdf/print.js`, through Playwright if it is installed, otherwise the Chromium command
line). Set `CHROME=/path/to/chrome` if Chromium is not found.

## How to cite

Wang, Q. (2026). *A wider window, the same wall: Two decades of item-level DIF tests for binary items*
(Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24447165

Later versions are listed on the Zenodo record. Code: https://github.com/wangqiu-syr/dif-power-binary-items

Citation metadata: [`CITATION.cff`](CITATION.cff).

## License

- Code: MIT ([`LICENSE`](LICENSE)).
- Report, supplement, figures and documentation: CC BY 4.0 ([`LICENSE-report`](LICENSE-report)).

## AI assistance

Code, analyses, documentation and drafts of the report were developed with Claude (Anthropic) as a
coding and writing assistant, under the author's direction. The author designed the study, checked
the results and is responsible for the content.

## References

Siepe, B. S., Bartoš, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W., & Pawel, S. (2024).
Simulation studies for methodological research in psychology: A standardized template for
planning, preregistration, and reporting. *Psychological Methods.* Advance online publication. https://doi.org/10.1037/met0000695

Vilhuber, L., Connolly, M., Koren, M., Llull, J., & Morrow, P. (2022). *A template README for
social science replication packages* (Version 1.1). Zenodo. https://doi.org/10.5281/zenodo.7293838
