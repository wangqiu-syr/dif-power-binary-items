## =============================================================================
## study_methods.R
## Methods added for the full power study (R/08_power_study.R; docs/power_study.md). Every item
## test returns a signed z (positive = easier for the focal group) and a two-sided p-value; for
## likelihood-ratio tests z = sign(estimate) * sqrt(LR).
## Uses R/power/pilot_methods.R (MH, cml_lr_one, the tree) and R/power/lrt_variants.R (difLRT).
## =============================================================================

signed <- function(z) list(z = z, p = 2 * pnorm(-abs(z)))

## CML Rasch fit. psychotools 0.7-3's analytic Hessian fails when a missing-data pattern has fewer
## than 4 items ("n > 1L is not TRUE" in elementary_symmetric_functions), e.g. a purified anchor of
## 2 items plus the split item. Only then refit with the numeric Hessian: the same estimates and
## log-likelihood; SEs agree with the analytic ones to 5 decimals where both run (decision D22).
rasch_cml <- function(r) tryCatch(psychotools::raschmodel(r), error = function(e)
  if (grepl("n > 1L", conditionMessage(e), fixed = TRUE)) psychotools::raschmodel(r, deriv = "numeric") else stop(e))

## ---- CML likelihood-ratio test, signed ------------------------------------------------------
## As cml_lr_one(), plus the sign of the estimated DIF (reference minus focal difficulty).
cml_lr_signed <- function(resp, x, use, j) {
  keep <- sort(union(use, j))
  ll0 <- as.numeric(logLik(rasch_cml(resp[, keep, drop = FALSE])))
  r2 <- cbind(resp[, setdiff(keep, j), drop = FALSE],
              jR = ifelse(x == 0, resp[, j], NA), jF = ifelse(x == 1, resp[, j], NA))
  m1 <- rasch_cml(r2)
  ip <- psychotools::itempar(m1)
  lr <- max(0, 2 * (as.numeric(logLik(m1)) - ll0))
  sign(ip[["jR"]] - ip[["jF"]]) * sqrt(lr)
}
cml_test <- function(resp, x, anchor_of) {
  signed(vapply(seq_len(ncol(resp)), function(j)
    tryCatch(cml_lr_signed(resp, x, anchor_of(j), j), error = function(e) NA_real_), 0))
}
## Standard purification: re-test every item against the unflagged items until the flagged set repeats
cml_purified_test <- function(resp, x, alpha = .05, nr_iter = 10) {
  all <- seq_len(ncol(resp)); r <- cml_test(resp, x, function(j) setdiff(all, j)); seen <- list(); it <- 0
  repeat {
    fl <- which(r$p < alpha)
    if (!length(fl) || it >= nr_iter || any(vapply(seen, identical, TRUE, fl))) break
    seen[[length(seen) + 1]] <- fl; it <- it + 1
    anc <- setdiff(all, fl); if (length(anc) < 2) break
    r <- cml_test(resp, x, function(j) setdiff(anc, j))
  }
  r
}

## ---- Robust anchor: the items closest to the median group difference ---------------------
## Fit the Rasch model by CML in each group separately; d_i = b_iR - b_iF (each group's item
## parameters centred at 0). Without DIF every d_i estimates the same constant (the impact shows
## up in the centring, not here); a DIF item is shifted by its DIF. The median of the d_i is not
## moved by DIF in fewer than half the items, so the k items closest to it make the anchor
## (excluding the tested item). Then the CML LRT against that anchor.
robust_anchor_d <- function(resp, x) {
  bR <- psychotools::itempar(psychotools::raschmodel(resp[x == 0, , drop = FALSE]))
  bF <- psychotools::itempar(psychotools::raschmodel(resp[x == 1, , drop = FALSE]))
  unname(bR - bF)
}
cml_robust_test <- function(resp, x, k = 5) {
  d <- tryCatch(robust_anchor_d(resp, x), error = function(e) NULL)
  all <- seq_len(ncol(resp))
  if (is.null(d) || anyNA(d) || any(!is.finite(d)))            # fall back: all other items
    return(c(cml_test(resp, x, function(j) setdiff(all, j)), fallback = TRUE))
  dist <- abs(d - median(d))
  anchor_of <- function(j) { o <- setdiff(order(dist), j); sort(o[seq_len(k)]) }
  c(cml_test(resp, x, anchor_of), fallback = FALSE)
}

## ---- Bundle test: one common shift on pre-specified items -----------------------------------
## The suspected items D are split into reference and focal virtual items; the other items are
## the anchor. delta_j = b_jR - b_jF with covariance V from the CML fit; the GLS common shift
## has z = 1'V^-1 delta / sqrt(1'V^-1 1).
bundle_test <- function(resp, x, D) {
  A <- setdiff(seq_len(ncol(resp)), D)
  vi <- do.call(cbind, lapply(D, function(j) cbind(ifelse(x == 0, resp[, j], NA), ifelse(x == 1, resp[, j], NA))))
  colnames(vi) <- as.vector(rbind(paste0("R", D), paste0("F", D)))
  r2 <- cbind(resp[, A, drop = FALSE], vi)
  ip <- psychotools::itempar(psychotools::raschmodel(r2)); V <- vcov(ip)
  C <- matrix(0, length(D), length(ip), dimnames = list(NULL, names(ip)))
  for (i in seq_along(D)) { C[i, paste0("R", D[i])] <- 1; C[i, paste0("F", D[i])] <- -1 }
  dl <- drop(C %*% ip); Vd <- C %*% V %*% t(C); w <- solve(Vd, rep(1, length(D)))
  z <- sum(w * dl) / sqrt(sum(w))
  list(z = z, p = 2 * pnorm(-abs(z)), shift = sum(w * dl) / sum(w))
}

## ---- MIMIC ------------------------------------------------------------------------------------
## test = "none": the chi-square test is not needed for the z of a direct effect (identical z, ~35% faster)
mimic_fit_fast <- function(model, dat, items) {
  lavaan::sem(model, data = dat, ordered = items, estimator = "WLSMV", parameterization = "delta",
              std.lv = TRUE, test = "none")
}
## The 2011 model: all ten direct effects free, no f ~ x
mimic_2011 <- function(dat, items) {
  m <- paste(paste("f1 =~", paste(items, collapse = " + ")), paste(paste0(items, " ~ x"), collapse = "\n"), sep = "\n")
  f <- tryCatch(suppressWarnings(mimic_fit_fast(m, dat, items)), error = function(e) NULL)
  if (is.null(f)) return(signed(rep(NA_real_, length(items))))
  signed(unname(mimic_z(f, items)))
}
## MIMIC with scale purification (M-SP; Wang, Shih & Yang, 2009): f ~ x; test each item with the
## unflagged items as the anchor; free the flagged items' direct effects; repeat until the flagged
## set repeats (at most max_iter passes).
mimic_msp <- function(dat, items, alpha = .05, max_iter = 4) {
  C <- length(items); meas <- paste("f1 =~", paste(items, collapse = " + "))
  fit_z <- function(free, which) {
    m <- paste(meas, "f1 ~ x", paste(paste0(items[free], " ~ x"), collapse = "\n"), sep = "\n")
    f <- tryCatch(suppressWarnings(mimic_fit_fast(m, dat, items)), error = function(e) NULL)
    if (is.null(f)) rep(NA_real_, length(which)) else unname(mimic_z(f, items[which]))
  }
  flagged <- integer(0); seen <- list(); it <- 0; z <- rep(NA_real_, C)
  repeat {
    it <- it + 1; z <- rep(NA_real_, C)
    for (j in setdiff(seq_len(C), flagged)) z[j] <- fit_z(sort(c(flagged, j)), j)
    if (length(flagged)) z[flagged] <- fit_z(flagged, flagged)
    new <- which(2 * pnorm(-abs(z)) < alpha)
    if (identical(new, flagged) || it >= max_iter || length(new) >= C - 1 ||
        any(vapply(seen, identical, TRUE, new))) break
    seen[[length(seen) + 1]] <- new; flagged <- new
  }
  c(signed(z), passes = it)
}
