# Supplement A3. Retrospective ADEMP-PreReg form for Simulation 2

[Qiu Wang](https://orcid.org/0000-0003-1635-8050) (Syracuse University) · Supplement to *A wider window, the same wall: Two decades of item-level DIF tests for binary items* · Version 1.0.0 · November 2026

Text: CC BY 4.0. Repository: https://github.com/wangqiu-syr/dif-power-binary-items

## About this form

This form follows the structure of the ADEMP-PreReg template for simulation studies (Siepe et
al., 2024). It was completed **retrospectively**, in November 2026, after the study had run. It is
**not a preregistration**. The plan was fixed in version control, in three commits of the author's
development repository, and not in a public registry. That repository is not public; the commits
can be shown to editors and reviewers on request, and the texts that matter are reproduced below
verbatim (§7).

Each answer gives the plan, the commit that fixed it, and what was done. Every difference between
the plan and the study is listed in §6.4.

**Table A3.1**

*The Commits That Fixed the Plan*

| Commit | Time (UTC) | What it fixed |
|---|---|---|
| `e2ae5a1` | 2026-11-19 00:32 | Factors and levels, methods, outcomes, 500 replications, H1–H5 (`docs/power.md` §7). Before any Simulation 2 data were generated |
| `4b26331` | 2026-11-19 01:06 | The design as run, the methods, the flag rules, H1–H6 (`docs/power_study.md` §2–4) and the code (`R/08_power_study.R`, `R/09_power_study_summary.R`, `R/power/study_methods.R`), committed together with the 20-replication pilot's results. Before the full run |
| `29dc42a` | 2026-11-19 09:33 | The results of the full run (40 conditions × 300 replications) and the decisions on H1–H6 |

*Note.* Full hashes: `e2ae5a1dfc4a89f79aa14692328eccf9e095a2f2`,
`4b263318ad8b65f71fdb0b931db83490d4d267bc`, `29dc42ae72cd45a2540e5da75b59c8df49ff2dd6`. Two commits fall between `4b26331` and `29dc42a`: `00dbc4b` (03:09 UTC, the robust-anchor
diagnostic `R/10_robust_anchor_check.R`) and `57ad383` (05:09 UTC, the layout of the result
figures). Neither changed the design or the hypotheses: §2–4 of `docs/power_study.md` are
identical at `4b26331` and `29dc42a`.

**Limits of this evidence.** Git commit times are set by the committing machine's clock. GitHub
kept a server-side record of the push of `29dc42a` (2026-11-19 09:33:23 UTC) but not of the two
earlier commits. The hash chain fixes the order of the commits; the clock times are
self-documented. In the terms of the TOP guidelines (Nosek et al., 2015) this is disclosure, not registration.

## 0. General information

**0.1 Title.** Simulation 2 of *A wider window, the same wall: Two decades of item-level DIF tests for binary items*: the
pre-specified power study.

**0.2 Contributors.** Qiu Wang (Syracuse University) designed the study, checked the results and
is responsible for the content. Code and drafts were developed with Claude (Anthropic) as a coding
and writing assistant.

**0.3 Description.** Simulation 1 and exact calculations showed that item-level DIF tests of 10
Rasch items with 50 examinees per group have power near .35 at best, and that the three methods
studied lose power or Type I error control in different ways. Simulation 2 tests the proposed
remedies where they matter: larger samples, DIF in both directions, groups that differ in ability
(impact), and a smaller DIF.

**0.4 Did the author develop any of the methods studied?** No published method. The
robust-anchor CML likelihood-ratio test (the k items closest to the median group difference) is a
simple version of the robust-scaling idea (Halpin, 2024), set up for this study and not tuned on
its results (decision D18 in `docs/decisions.md`).

## 1. Aims

**1.1 Aim.** Plan (`e2ae5a1`): "A design to test the levers where they matter." As run
(`4b26331`): test the remedies for low power and for Type I error inflation identified in Section 4
of the report (sample size, a valid anchor, a correctly computed LRT, a one-sided test, a joint
test of suspected items, dropping the tree's gate, false discovery rate control) against six
hypotheses fixed before the results (§5.5).

## 2. Data-generating mechanism

**2.1 How are the parameters specified?** Parametric. The Rasch model with the 10 difficulties of
the baseline design, $b = (-1, -1, -1, 0, 0, 0, 0, 1, 1, 1)$. Uniform DIF on i1, i5 and i10. θ is
normal in each group: $\theta_R \sim N(0, 1)$; $\theta_F \sim N(0, 1)$ or $N(-0.5, 1)$. Responses are
drawn as in Simulation 1 (`respon_1pl()`, in the draw order of the original baseline study). Planned at `e2ae5a1`; specified in full at
`4b26331`.

**2.2 Factors.** Examinees per group; DIF direction; impact; DIF size. Unchanged from `e2ae5a1`.

**2.3 Levels.**
- Examinees per group: 50, 100, 150, 250 (equal groups).
- DIF direction: one-way (i1, i5, i10 all easier for the focal group); mixed (i1 and i10 easier,
  i5 harder).
- Impact: none; or the focal group 0.5 SD lower.
- DIF size: 0.5 or 0.8 logits.

Unchanged from `e2ae5a1`.

**2.4 How are the levels combined?** Full factorial for the DIF conditions, 4 × 2 × 2 × 2 = 32,
plus 8 no-DIF conditions (each sample size × impact), for 40 conditions. Unchanged from `e2ae5a1`.

## 3. Estimands and other targets

**3.1 Targets.** For each item j, the null hypothesis of no DIF ($\delta_j = 0$). The methods are
evaluated as tests. Plan (`e2ae5a1`) also listed the bias of the DIF estimates as an outcome; it
was not estimated (§6.4, deviation 4).

## 4. Methods

**4.1 Which methods, and which quantities are extracted?**

**Table A3.2**

*Methods: Plan and As Run*

| Method | Plan (`e2ae5a1`) | As run (`4b26331`) |
|---|---|---|
| Benchmark: θ known (logistic regression with θ as an offset) | not listed | run |
| Benchmark: CML LRT with the true clean anchor | not listed | run |
| MIMIC, baseline specification | reference | run |
| MIMIC with `f1 ON x` and scale purification (M-SP; Wang et al., 2009) | 4 | run (at most 4 passes) |
| Mantel–Haenszel, purified, no continuity correction | 3 | run |
| CML LRT, standard purification | 2 | run |
| CML LRT, robust anchor | 1 | run with k = 5 and k = 3 |
| Rasch tree (gate) + anchortest | reference | run |
| anchortest without the gate | 5 | run |
| difLRT as in the baseline study | reference, "at 50 per group only" | run in two n = 50, no-impact conditions, one-way DIF of 0.8 and no DIF (the lme4 rebuild of difR 6.1.0's statistic) |
| Bundle test of i1 + i5 + i10 | 6: against the purified anchor, one-way conditions only | the other seven items as the anchor, every condition |
| One-sided α = .05 and Benjamini–Hochberg across the 10 items | 7 | applied afterwards to the saved statistics |

Quantities saved for every method, item and replication: the signed statistic ($z$, or $\mathrm{sign} \times \sqrt{LR}$)
and the two-sided p-value; for the bundle, its z and p. Per replication: the tree's global p-value,
the number of M-SP passes, whether the robust anchor fell back to all other items, and the run
time. All are in `results/power_study/replications.rds` and `replications_info.csv`.

## 5. Performance measures

**5.1 Which measures?** Plan (`e2ae5a1`): "Power and Type I error per method and cell;
size-adjusted power; the false discovery rate; the bias of the DIF estimates." As run: power (share
of tests on i1, i5 and i10 that flag) and Type I error (share of tests on the clean items that flag:
seven items, or all ten in the no-DIF conditions), at two-sided α = .05 (primary); one-sided power
and Type I error in the hypothesized direction (one-way and no-DIF conditions); Benjamini–Hochberg
power and false discovery rate; size-adjusted power (critical p-value from the no-DIF condition with
the same n and impact); fit failures. Not done: bias.

**5.2 How is Monte Carlo uncertainty computed and reported?** Not specified in the plan, beyond the
precision targets in 5.3. The report gives, for every rate, an MCSE that treats the replication as
the cluster (items in one data set are correlated):
$\text{MCSE}(R) = \sqrt{\mathrm{Var}(f_r - R\, t_r) / n_{\text{rep}}} \,/\, \bar t$,
where $R$ is the pooled rate, replication $r$ has $f_r$ flags among $t_r$ tests, and $\bar t$ is
the mean number of tests per replication. Computed after the study
by `R/report_tables.R` (report, Appendix A2).

**5.3 How many replications?** Plan (`e2ae5a1`): 500 per condition ("power to ±.02, Type I to
±.01"). Changed at `4b26331`, before the full run: 300 per condition, because the pilot's timing
(about 4 s per replication at 50 per group and 12 s at 250) put 500 at about 13 hours on 4 cores
(decision D17). The pilot's 20 replications per condition are replications 1–20 of the 300.
Achieved precision: largest MCSE .024 for power and .014 for Type I error.

**5.4 How are non-convergence and missing values handled?** A failed fit returns NA and is left out
of the denominator, as in Simulation 1 (decision D10). If the group-wise CML fits behind the robust
anchor fail, it falls back to all other items as the anchor, and the fallback is counted. The first
anchor item of anchortest fixes the scale and is not tested; it is left out of the denominator.
Specified in the code at `4b26331`. As run: no fit failed in any method and the robust anchor never
fell back (`results/power_study/fit_problems.csv`).

**5.5 How are the results interpreted?** Against six hypotheses. H1–H5 were fixed at `e2ae5a1`; at
`4b26331`, H2, H4 and H5 were extended and H6 was added (exact texts in §7). The decisions, recorded
with the results at `29dc42a`:

**Table A3.3**

*Decisions on the Hypotheses*

| Hypothesis (as at `4b26331`) | Decision | Evidence |
|---|---|---|
| H1. Robust-anchor CML LRT: Type I within .05 ± .015 in every cell; power within .03 of the clean-anchor ceiling | Rejected | Power on average .077 below the ceiling (criterion .03); Type I .066–.080 in 8 of 16 one-way conditions, the smaller exceedances within Monte Carlo error |
| H2. One-way DIF inflates the tests whose first-pass anchor is all other items, growing with n; purification removes most of the excess only at n ≥ 150; mixed DIF cancels | Partly supported | Inflation .05–.09 under one-way DIF; mixed DIF cancels (.042–.064); the prediction "only at n ≥ 150" was wrong: purification removes most of the excess at every n |
| H3. Baseline MIMIC: nominal without impact, near 1 with impact | Partly supported | .046–.061 without impact; .17 / .29 / .42 / .60 with impact (the two-proportion prediction); the prediction "near 1" was wrong |
| H4. Bundle test: .80 at about 64 per group (0.8 logits, one-way); holds .05 without DIF | Supported | .70 at 50 and .95 at 100; .027–.077 without DIF (mean .049) |
| H5. Valid anchor: .80 at 140–180 per group for 0.8 logits; below .80 at 250 for 0.5 logits | Supported | Clean-anchor benchmark about 150–160 per group; .61–.62 at 250 for 0.5 logits |
| H6. The gate costs power at every n; ungated anchortest above .05 under one-way DIF | Supported | Gate cost .14 / .12 / .09 / .05 (mean power over the DIF conditions); ungated Type I .073–.117 |

## 6. Other

**6.1 Software.** R 4.3.3; lavaan 0.6-17 (MIMIC, WLSMV); psychotools 0.7-3 (CML Rasch fits,
anchortest); psychotree 0.16-0 and partykit 1.2-20 (Rasch tree); lme4 1.1-35.1 (the difLRT
rebuild); base R (Mantel–Haenszel, logistic regression). Versions as run:
`results/power_study/sessionInfo.txt`.

**6.2 Computational environment.** Linux (Ubuntu 24.04), 4 cores. The cloud machine running the
study was reclaimed twice, so the run finished in three legs. Each replication has its own seed and
its own output file, and a restart skips finished files, so the results are those of an
uninterrupted run.

**6.3 Reproducibility.** Replication r of condition c uses `set.seed(20261119 + 1000c + r)`;
anchortest's tie-breaking uses that seed + 7. Every statistic is saved, so other flag rules can be
applied without rerunning. `tests/run_tests.R` regenerates replication 1 of condition 1 and
compares it with the stored statistics; `R/run_all.R FULL` reruns the whole study.

**6.4 Deviations and other disclosures.**

**Table A3.4**

*Deviations From the Plan*

| # | Plan | As run | When decided |
|---|---|---|---|
| 1 | 500 replications per condition | 300 | `4b26331`, from the pilot's timing, before the full run |
| 2 | difLRT "at 50 per group only" | Only the two n = 50, no-impact conditions (one-way DIF of 0.8 and no DIF), where it is the baseline reference | `4b26331`, before the full run |
| 3 | Bundle test against the purified anchor, one-way conditions only | The other seven items as the anchor; every condition | `4b26331`, before the full run |
| 4 | Bias of the DIF estimates as an outcome | Not estimated: the estimates were not saved. Part of the follow-up study (https://osf.io/35nvh/) | at the analysis |
| 5 | H2, H4, H5 as at `e2ae5a1` | Extended at `4b26331` (texts in §7) | `4b26331` |
| 6 | No H6 | H6 added | `4b26331` |
| 7 | Two benchmarks not listed | θ known and the CML LRT with the true clean anchor added | `4b26331` |

Other disclosures:
- **The pilot.** Deviations 1–3 and 5–7 were committed together with the results of the
  20-replication pilot (`4b26331`), so they may have been informed by it. The pilot's data sets are
  replications 1–20 of the full run, not an independent sample.
- **During the run.** `R/10_robust_anchor_check.R` (commit `00dbc4b`, 03:09 UTC) was written while
  the full run was in progress, after results from the first conditions to finish were available
  (the run scheduled the two difLRT conditions and then n = 250 first; the check examines
  condition 25, n = 250). It explains a result and tests nothing; the report labels it exploratory.
  The figure layout was changed during the run as well (`57ad383`).
- **Analyses added after the results.** The replication-level MCSEs (§5.2) and the report's
  figures and tables were produced after `29dc42a` from the saved results. They change no rate.

## 7. The hypotheses as committed

In the quoted texts, "2011" refers to the baseline study and its specifications (the report's
Sections 1–3); the texts are left as committed.

**At `e2ae5a1`** (`docs/power.md` §7, "Hypotheses, from sections 2–3 and 5"), verbatim:

> - H1. The CML LRT with a robust anchor holds Type I error within .05 ± .015 in every cell, and its
>   power is within .03 of the clean-anchor ceiling.
> - H2. With one-way DIF, the all-other-anchor tests' Type I error grows with n (predicted .08, .12,
>   .15, .22 at 50, 100, 150, 250 for 0.8 logits); purification removes most of it only at n ≥ 150.
> - H3. The 2011 MIMIC's Type I error is nominal without impact and near 1 with impact.
> - H4. The bundle test reaches .80 at about 64 per group for 0.8 logits.
> - H5. Power at 0.8 logits with a valid anchor reaches .80 at 140–180 per group (section 2).

**At `4b26331`** (`docs/power_study.md` §4, "Hypotheses (from power.md, fixed before the
results)"), verbatim:

> - **H1.** The robust-anchor CML LRT holds Type I error within .05 ± .015 in every cell, and its
>   power is within .03 of the clean-anchor ceiling.
> - **H2.** With one-way DIF, tests that use all other items as the anchor in their first pass
>   (purified CML LRT, purified MH, M-SP) have Type I error above .05 that grows with n (the
>   exact calculation for the unpurified test: .08, .12, .15, .22 at 50, 100, 150, 250 for 0.8
>   logits). Purification removes most of the excess only at n ≥ 150. With mixed DIF the
>   contamination largely cancels.
> - **H3.** The 2011 MIMIC's Type I error is nominal without impact and near 1 with impact.
> - **H4.** The bundle test reaches .80 at about 64 per group for 0.8 logits in the one-way cells,
>   and holds .05 in the no-DIF cells.
> - **H5.** With a valid anchor, power for 0.8 logits reaches .80 at 140–180 per group, and for 0.5
>   logits stays below .80 even at 250 (it needs about 350).
> - **H6.** The Rasch tree's gate costs power at every n, and the ungated anchortest's Type I error is
>   above .05 under one-way DIF.

**What changed between the two commits:** H2 names the tests, says "above .05", and adds the
mixed-DIF prediction; H4 adds "in the one-way cells" and "holds .05 in the no-DIF cells"; H5 adds
the 0.5-logit prediction; H6 is new. H1 and H3 are unchanged.

## References

Halpin, P. F. (2024). Differential item functioning via robust scaling. *Psychometrika, 89*(3), 796–821. https://doi.org/10.1007/s11336-024-09957-6

Nosek, B. A., Alter, G., Banks, G. C., Borsboom, D., Bowman, S. D., Breckler, S. J., Buck, S., Chambers, C. D., Chin, G., Christensen, G., Contestabile, M., Dafoe, A., Eich, E., Freese, J., Glennerster, R., Goroff, D., Green, D. P., Hesse, B., Humphreys, M., . . . Yarkoni, T. (2015). Promoting an open research culture. *Science, 348*(6242), 1422–1425. https://doi.org/10.1126/science.aab2374

Siepe, B. S., Bartoš, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W., & Pawel, S. (2024). Simulation studies for methodological research in psychology: A standardized template for planning, preregistration, and reporting. *Psychological Methods.* Advance online publication. https://doi.org/10.1037/met0000695

Wang, W.-C., Shih, C.-L., & Yang, C.-C. (2009). The MIMIC method with scale purification for detecting differential item functioning. *Educational and Psychological Measurement, 69*(5), 713–731. https://doi.org/10.1177/0013164409332228
