## =============================================================================
## pilot_methods.R
## Methods for the power pilot (R/04_power_pilot.R; docs/power.md, section 5). Each returns, per
## item, a statistic and a two-sided p-value (signed z where the method has a sign: positive =
## easier for the focal group, the direction of the simulated DIF).
## =============================================================================

## Data generation identical, draw for draw, to simulate_dif_data(), but also returning theta.
simulate_dif_theta <- function(b_ref, b_foc, n_ref = 50, n_foc = 50) {
  th_r <- rnorm(n_ref); ref <- respon_1pl(b_ref, th_r)
  th_f <- rnorm(n_foc); foc <- respon_1pl(b_foc, th_f)
  dat <- data.frame(c(rep(0L, n_ref), rep(1L, n_foc)), rbind(ref, foc))
  names(dat) <- c("x", paste0("i", seq_along(b_ref)))
  list(dat = dat, theta = c(th_r, th_f))
}

z_out <- function(z) list(z = z, p = 2 * pnorm(-abs(z)))

## ---- Benchmarks: what the data could give at best ------------------------------------------
## Oracle: theta known, Rasch slope known (= 1): logit P = theta - b_j + delta_j x.
oracle_theta <- function(resp, x, theta) {
  z_out(vapply(seq_len(ncol(resp)), function(j) {
    f <- glm(resp[, j] ~ x, family = binomial, offset = theta)
    unname(coef(summary(f))["x", "z value"]) }, 0))
}
## Marginal: compare the groups' proportions correct, assuming no impact (logit Wald test).
marginal_test <- function(resp, x) {
  z_out(vapply(seq_len(ncol(resp)), function(j) {
    f <- glm(resp[, j] ~ x, family = binomial)
    unname(coef(summary(f))["x", "z value"]) }, 0))
}

## ---- Mantel-Haenszel (the Rasch-efficient conditional test; Holland & Thayer, 1988) --------
## Matching score = sum of the anchor items + the studied item. Continuity correction as in
## difR::difMH(). Signed by the MH common odds ratio (positive z = focal favoured).
mh_one <- function(yj, x, score, correct = TRUE) {
  A <- E <- V <- num <- den <- 0
  for (s in unique(score)) {
    k <- score == s; Tk <- sum(k); if (Tk < 2) next
    nR <- sum(k & x == 0); nF <- Tk - nR; m1 <- sum(yj[k]); m0 <- Tk - m1
    if (nR == 0 || nF == 0 || m1 == 0 || m0 == 0) next
    a <- sum(yj[k & x == 0]); b <- nR - a; cc <- sum(yj[k & x == 1]); d <- nF - cc
    A <- A + a; E <- E + nR * m1 / Tk; V <- V + as.numeric(nR) * nF * m1 * m0 / (Tk^2 * (Tk - 1))   # double: integer counts overflow above ~430 per score level
    num <- num + a * d / Tk; den <- den + b * cc / Tk
  }
  if (V == 0) return(NA_real_)
  chi <- (max(abs(A - E) - if (correct) .5 else 0, 0))^2 / V
  -sign(log(num / den)) * sqrt(chi)          # log OR > 0: reference favoured
}
mh_test <- function(resp, x, anchor_of = function(j) setdiff(seq_len(ncol(resp)), j), correct = TRUE) {
  z_out(vapply(seq_len(ncol(resp)), function(j) {
    mh_one(resp[, j], x, rowSums(resp[, c(anchor_of(j), j), drop = FALSE]), correct) }, 0))
}
## Purification as usually defined: re-test every item against the unflagged items.
mh_purified <- function(resp, x, alpha = .05, nr_iter = 10, correct = TRUE) {
  C <- ncol(resp); r <- mh_test(resp, x, correct = correct); seen <- list(); it <- 0
  repeat {
    fl <- which(r$p < alpha)
    if (!length(fl) || it >= nr_iter || any(vapply(seen, identical, TRUE, fl))) break
    seen[[length(seen) + 1]] <- fl; it <- it + 1
    anc <- setdiff(seq_len(C), fl)
    r <- mh_test(resp, x, function(j) setdiff(anc, j), correct)
  }
  r
}

## ---- MIMIC variants (lavaan, WLSMV, delta parameterization, factor variance 1) ---------------
mimic_fit <- function(model, dat, items) {
  lavaan::sem(model, data = dat, ordered = items, estimator = "WLSMV", parameterization = "delta",
              std.lv = TRUE)
}
mimic_z <- function(fit, which) {
  pe <- lavaan::parameterEstimates(fit); pe <- pe[pe$op == "~" & pe$rhs == "x", ]
  setNames(pe$z, pe$lhs)[which]
}
## (a) the 2011 model: all ten direct effects free, no f ~ x (equal group means assumed)
## (b) the same with equal loadings (the generating model is Rasch)
## (c) anchored: f ~ x, only the tested item's direct effect free (all other items are the anchor)
## (d) oracle anchor: f ~ x, direct effects free for every item outside the true clean set
mimic_variants <- function(dat, items, clean) {
  C <- length(items); all <- seq_len(C)
  meas <- paste("f1 =~", paste(items, collapse = " + "))
  meas_eq <- paste("f1 =~", paste0("l*", items, collapse = " + "))
  de <- function(js) paste(paste0(items[js], " ~ x"), collapse = "\n")
  safe <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) NULL)
  out <- list()
  f <- safe(mimic_fit(paste(meas, de(all), sep = "\n"), dat, items))
  out$mimic2011 <- z_out(if (is.null(f)) rep(NA, C) else unname(mimic_z(f, items)))
  f <- safe(mimic_fit(paste(meas_eq, de(all), sep = "\n"), dat, items))
  out$mimic_eqload <- z_out(if (is.null(f)) rep(NA, C) else unname(mimic_z(f, items)))
  out$mimic_anchored <- z_out(vapply(all, function(j) {
    f <- safe(mimic_fit(paste(meas, "f1 ~ x", de(j), sep = "\n"), dat, items))
    if (is.null(f)) NA_real_ else unname(mimic_z(f, items[j])) }, 0))
  dirty <- setdiff(all, clean)
  zo <- rep(NA_real_, C)
  if (length(dirty)) {
    f <- safe(mimic_fit(paste(meas, "f1 ~ x", de(dirty), sep = "\n"), dat, items))
    if (!is.null(f)) zo[dirty] <- mimic_z(f, items[dirty])
  }
  for (j in clean) {
    f <- safe(mimic_fit(paste(meas, "f1 ~ x", de(sort(c(dirty, j))), sep = "\n"), dat, items))
    if (!is.null(f)) zo[j] <- mimic_z(f, items[j])
  }
  out$mimic_oracle <- z_out(zo)
  out
}

## ---- Rasch tree: the global test's p-value, and anchortest() with and without the gate --------
tree_variants <- function(resp, x, items, tie_seed) {
  if (!"package:psychotree" %in% search()) suppressPackageStartupMessages(library(psychotree))
  C <- ncol(resp); dd <- data.frame(x = factor(x)); dd$resp <- resp
  na <- list(z = rep(NA_real_, C), p = rep(NA_real_, C))
  root_p <- tryCatch({
    tr <- psychotree::raschtree(resp ~ x, data = dd)
    partykit::info_node(partykit::node_party(tr))$test["p.value", 1]
  }, error = function(e) NA_real_)
  set.seed(tie_seed)
  at <- tryCatch(suppressWarnings(psychotools::anchortest(resp ~ x, data = dd, class = "forward",
                                                          select = "MTT")), error = function(e) NULL)
  ungated <- na
  if (!is.null(at)) {
    tested <- as.integer(sub("^.*?(\\d+)$", "\\1", names(at$final_tests$test$tstat)))
    ungated$z[tested] <- at$final_tests$test$tstat; ungated$p[tested] <- at$final_tests$test$pvalues
  }
  list(tree_root_p = root_p, anchortest = ungated)
}

## ---- Conditional-ML likelihood-ratio test (a correctly computed Rasch LRT) -------------------
## Compact: the Rasch model on the items in `use` plus item j. Augmented: item j split into a
## reference and a focal "virtual item" (psychotools::raschmodel handles the structural NAs).
## Conditioning on the raw score removes theta and its distribution, so impact needs no model.
cml_lr_one <- function(resp, x, use, j) {
  keep <- sort(union(use, j))
  ll0 <- as.numeric(logLik(psychotools::raschmodel(resp[, keep, drop = FALSE])))
  r2 <- cbind(resp[, setdiff(keep, j), drop = FALSE],
              jR = ifelse(x == 0, resp[, j], NA), jF = ifelse(x == 1, resp[, j], NA))
  max(0, 2 * (as.numeric(logLik(psychotools::raschmodel(r2))) - ll0))
}
cml_lr <- function(resp, x, anchor_of) {
  st <- vapply(seq_len(ncol(resp)), function(j)
    tryCatch(cml_lr_one(resp, x, anchor_of(j), j), error = function(e) NA_real_), 0)
  list(z = st, p = pchisq(st, 1, lower.tail = FALSE))          # z holds the LR statistic
}
cml_variants <- function(resp, x, clean, alpha = .05, nr_iter = 10) {
  C <- ncol(resp); all <- seq_len(C)
  none <- cml_lr(resp, x, function(j) setdiff(all, j))
  r <- none; seen <- list(); it <- 0
  repeat {                                                     # standard purification
    fl <- which(r$p < alpha)
    if (!length(fl) || it >= nr_iter || any(vapply(seen, identical, TRUE, fl))) break
    seen[[length(seen) + 1]] <- fl; it <- it + 1
    anc <- setdiff(all, fl); if (length(anc) < 2) break
    r <- cml_lr(resp, x, function(j) setdiff(anc, j))
  }
  list(cml_none = none, cml_purified = r, cml_oracle = cml_lr(resp, x, function(j) setdiff(clean, j)))
}
