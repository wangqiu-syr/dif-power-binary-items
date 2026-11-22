# difLRT: the IRT likelihood-ratio DIF test

## Its role in this study

`difLRT` is **not a cross-check**. The study compares three DIF detection methods on the same
simulated data sets:

1. **MIMIC**: a one-factor model with each item regressed on group (Mplus in 2011, lavaan here)
2. **IRT likelihood-ratio test**: `difLRT()` from the R package **difR**
3. **Rasch tree**: `raschtree()` from the R package **psychotree**, added in this rebuild
   ([methods.md §5](methods.md#5-rasch-tree-psychotree-added-in-this-rebuild))

The question is which method finds the true DIF items more often (power) while wrongly flagging
the clean items less often (Type I error). This note covers difLRT, and compares it with MIMIC,
the method it was paired with in 2011.

The cross-check in this replication is something else: the lavaan MIMIC results compared with
the 2011 Mplus 4.2 output (`R/01_single_dataset.R`). That shows lavaan is a faithful stand-in for Mplus.

## What difLRT does

`difLRT()` is from difR (Magis, Béland, Tuerlinckx, & De Boeck, 2010). It applies the
likelihood-ratio test for DIF (Thissen, Steinberg, & Wainer, 1988) under the **Rasch model**,
which gives each item a single difficulty parameter.

For each item *j* it fits two models to the reference and focal groups together:

| Model | Item *j*'s difficulty | Other (anchor) items |
|---|---|---|
| Compact | equal in both groups | equal in both groups |
| Augmented | **free to differ** between groups | equal in both groups |

The test statistic is the drop in deviance:

```
LRT_j = -2 * [ logL(compact) - logL(augmented) ]
```

If item *j* has no DIF, `LRT_j` follows a χ² distribution with 1 df. The item is flagged when
`LRT_j >= thr`, where `thr = qchisq(1 - alpha, 1)`, which is **3.84** at α = .05.

Because the Rasch model has no discrimination parameter, difLRT tests **uniform DIF** only: the
item is uniformly easier or harder for one group at every ability level. That is exactly what
this design simulates. Items i1, i5 and i10 are 0.8 logits easier for the focal group.

The anchor items are held equal across groups, which links the two groups' scales. So a real
difference in the groups' average ability (impact) is not mistaken for DIF.

## Purification (`purify = TRUE`)

If some anchor items have DIF themselves, the link between the groups is distorted, and that
can hide DIF or create false DIF. Purification reduces this:

1. Test every item, using all the other items as anchors.
2. Remove the flagged items from the anchor set and test every item again.
3. Repeat until the flagged set stops changing, or until `nrIter` iterations (default 10).

The 2011 code and this replication both use `purify = TRUE`.

## The call in this repository

From `run_lrt()` in `R/functions.R`, the same call as the 2011 code:

```r
res <- difR::difLRT(Data = data[, items],     # 0/1 item responses, i1-i10
                    group = data$x,           # 0 = reference, 1 = focal
                    focal.name = 1,
                    alpha = 0.05,
                    purify = TRUE)

res$LRT                 # LRT statistic for each item
res$thr                 # critical value, 3.84 for alpha = .05
res$LRT >= res$thr      # DIF flags (stored as LRTstat in 02_simulation.R)
res$DIFitems            # which items were flagged
```

For the full argument list, the other output components, and how difR estimates the Rasch
model, see `?difLRT` in R.

## How the two methods differ in this study

| | difLRT | MIMIC (as specified in 2011) |
|---|---|---|
| Item model | Rasch: logistic, equal discrimination | one-factor probit: loadings free per item, like a 2PL |
| Matches the data-generating model? | yes, the data are 1PL | more general; estimates extra parameters |
| Estimation | maximum likelihood | WLSMV |
| DIF test | likelihood ratio, χ²(1) | Wald z for `item ON x` |
| 2011 flag rule | `LRT >= 3.84`: two-sided, α = .05 | `z >= 1.96`: one-sided, α = .025 |
| Group mean difference | allowed; the anchors fix the scale | fixed at 0 (no `f1 ON x`) |
| Anchor | purified, iterative | none; all ten direct effects are free |

Two consequences when reading the results:

- **The flag rules use different α levels.** At face value, MIMIC is held to a stricter
  false-positive rate. The simulation saves the raw z and LRT values, so both methods can be
  re-scored at the same α.
- **Different assumptions about impact.** MIMIC's equal-means assumption holds in this design,
  where both groups are N(0, 1). If impact is added, the MIMIC model needs `f1 ~ x` and an anchor
  item; difLRT needs no change.

## References

Finch, H. (2005). The MIMIC model as a method for detecting DIF: Comparison with Mantel-Haenszel,
SIBTEST, and the IRT likelihood ratio. *Applied Psychological Measurement, 29*(4), 278–295.

Magis, D., Béland, S., Tuerlinckx, F., & De Boeck, P. (2010). A general framework and an R package
for the detection of dichotomous differential item functioning. *Behavior Research Methods, 42*(3), 847–862.

Thissen, D., Steinberg, L., & Wainer, H. (1988). Use of item response theory in the study of group
differences in trace lines. In H. Wainer & H. I. Braun (Eds.), *Test validity* (pp. 147–169). Erlbaum.
