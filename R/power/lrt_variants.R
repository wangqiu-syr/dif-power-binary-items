## =============================================================================
## lrt_variants.R
## The IRT likelihood-ratio DIF test of difR::difLRT(), rebuilt with lme4 so that its parts can
## be varied one at a time (docs/power.md, section 3). difR 6.1.0's LRT() fits the Rasch model as a
## GLMM: item fixed effects, a random person intercept, and a group main effect (impact):
##   compact:   y ~ items - 1 + member + (1 | person)
##   augmented: compact + member:X_j          (X_j = 1 on item j's rows)
## and LRT_j = deviance(compact) - deviance(augmented). All other items are the anchor.
## =============================================================================

## Long format once per data set: one row per person x item.
lrt_long <- function(resp, group) {
  N <- nrow(resp); C <- ncol(resp)
  data.frame(y = as.vector(as.matrix(resp)), person = factor(rep(seq_len(N), C)),
             item = factor(rep(seq_len(C), each = N)), member = rep(group, C))
}

## LRT for each item in `test`, with the anchor = the items in `use` other than the tested one.
## `use` are the columns kept in the data (difR's purification drops flagged items from the data).
lrt_stats <- function(resp, group, use = seq_len(ncol(resp)), test = use) {
  d <- lrt_long(resp[, use, drop = FALSE], group)
  ## One optimizer (bobyqa) instead of lme4's default bobyqa + Nelder-Mead: the same model and
  ## Laplace likelihood as difR, about 4x faster. On the checked data sets the statistics agree with
  ## the default optimizer to within about 0.01 (docs/power.md, section 3.1).
  ctl <- lme4::glmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
  m0 <- lme4::glmer(y ~ item - 1 + member + (1 | person), data = d, family = binomial, control = ctl)
  dev0 <- stats::deviance(m0)
  vapply(test, function(j) {
    d$X <- as.numeric(d$item == which(use == j))
    m1 <- lme4::glmer(y ~ item - 1 + member + member:X + (1 | person), data = d, family = binomial,
                      control = ctl)
    dev0 - stats::deviance(m1)
  }, 0)
}

## The variants. thr = qchisq(.95, 1).
##  "difR"     : difR's purification, verbatim: test all items; drop the flagged ones from the data;
##               re-test only the remaining items with the remaining items; flags accumulate and the
##               first-round statistic of a flagged item is never updated.
##  "none"     : no purification (all other items as the anchor).
##  "standard" : purification as usually defined: re-test EVERY item against the current anchor
##               (the unflagged items); items can be unflagged; iterate until the flagged set repeats.
##  "oracle"   : the true DIF-free items as the anchor (the upper bound for an anchor-based LRT).
lrt_variant <- function(resp, group, variant = c("difR", "none", "standard", "oracle"),
                        clean = NULL, nr_iter = 10, thr = stats::qchisq(.95, 1)) {
  variant <- match.arg(variant); C <- ncol(resp); all <- seq_len(C)
  if (variant == "none") return(lrt_stats(resp, group))
  if (variant == "oracle") {
    ## clean items: tested with the other clean items as the anchor; DIF items: tested with the
    ## clean items as the anchor (the data hold the clean items plus the tested item).
    s <- numeric(C)
    s[clean] <- lrt_stats(resp, group, use = clean)
    for (j in setdiff(all, clean)) s[j] <- lrt_stats(resp, group, use = sort(c(clean, j)), test = j)
    return(s)
  }
  if (variant == "difR") {
    s <- lrt_stats(resp, group)
    if (max(s) <= thr) return(s)
    nodif <- all[s <= thr]; it <- 0
    while (it < nr_iter && length(nodif) >= 2) {
      it <- it + 1
      s2 <- lrt_stats(resp, group, use = nodif)
      s[nodif] <- s2
      if (max(s2) <= thr) break
      nodif <- nodif[s2 <= thr]
    }
    return(s)
  }
  ## standard purification
  s <- lrt_stats(resp, group); flagged <- all[s > thr]; seen <- list(flagged); it <- 0
  while (length(flagged) && it < nr_iter) {
    it <- it + 1
    anchor <- setdiff(all, flagged)
    if (length(anchor) < 1) break
    s <- numeric(C)
    s[anchor] <- if (length(anchor) >= 2) lrt_stats(resp, group, use = anchor) else 0
    for (j in flagged) s[j] <- lrt_stats(resp, group, use = sort(c(anchor, j)), test = j)
    new <- all[s > thr]
    if (any(vapply(seen, identical, TRUE, new))) break
    seen[[length(seen) + 1]] <- new; flagged <- new
  }
  s
}

## All four variants on one data set, sharing the first round (the "none" statistics).
lrt_all <- function(resp, group, clean, nr_iter = 10, thr = stats::qchisq(.95, 1)) {
  C <- ncol(resp); all <- seq_len(C)
  s1 <- lrt_stats(resp, group)
  ## difR
  sd <- s1
  if (max(s1) > thr) {
    nodif <- all[s1 <= thr]; it <- 0
    while (it < nr_iter && length(nodif) >= 2) {
      it <- it + 1; s2 <- lrt_stats(resp, group, use = nodif); sd[nodif] <- s2
      if (max(s2) <= thr) break
      nodif <- nodif[s2 <= thr]
    }
  }
  ## standard
  ss <- s1; flagged <- all[s1 > thr]; seen <- list(flagged); it <- 0
  while (length(flagged) && it < nr_iter) {
    it <- it + 1; anchor <- setdiff(all, flagged)
    if (length(anchor) < 2) break
    ss <- numeric(C)
    ss[anchor] <- lrt_stats(resp, group, use = anchor)
    for (j in flagged) ss[j] <- lrt_stats(resp, group, use = sort(c(anchor, j)), test = j)
    new <- all[ss > thr]
    if (any(vapply(seen, identical, TRUE, new))) break
    seen[[length(seen) + 1]] <- new; flagged <- new
  }
  ## oracle anchor
  so <- numeric(C)
  if (length(clean) == C) so <- s1 else {
    so[clean] <- lrt_stats(resp, group, use = clean)
    for (j in setdiff(all, clean)) so[j] <- lrt_stats(resp, group, use = sort(c(clean, j)), test = j)
  }
  list(none = s1, difR = sd, standard = ss, oracle = so)
}
