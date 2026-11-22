## =============================================================================
## functions.R
## R replication of the 2011 NCME DIF study: MIMIC (Mplus, WLSMV) vs. IRT
## likelihood-ratio test (difR::difLRT) for detecting uniform DIF, plus a third
## method added in this rebuild: the Rasch tree (psychotree).
##
## Mplus is replaced by lavaan (same model, same estimator). Data generation is
## the original 1PL code from original/Mplus_DIFR_R.r.
## =============================================================================

ITEMS <- paste0("i", 1:10)

## ---- 1. Data generation ------------------------------------------------------

## 1PL (Rasch) responses: n x k matrix of 0/1.
## Same model and same random-number draws, in the same order, as the 2011
## respon.1pl(): one runif(k * n) vector compared with P = 1 / (1 + exp(-(theta - b))),
## filled column by column (persons within items).
respon_1pl <- function(b, theta) {
  b <- as.vector(b)
  n <- length(theta)
  k <- length(b)
  p <- plogis(outer(theta, b, "-"))          # n x k, P(correct)
  matrix(as.integer(runif(k * n) < p), nrow = n, ncol = k)
}

## One replication's data set, laid out like the 2011 DIF_data:
## column x (0 = reference, 1 = focal) followed by items i1..ik.
## Reference group is generated first, then focal (the 2011 order), with
## theta ~ N(0, 1) in both groups (no impact).
simulate_dif_data <- function(b_ref, b_foc, n_ref = 50, n_foc = 50) {
  ref <- respon_1pl(b_ref, rnorm(n_ref, 0, 1))   # rgroup()
  foc <- respon_1pl(b_foc, rnorm(n_foc, 0, 1))   # fgroup()
  dat <- data.frame(c(rep(0L, n_ref), rep(1L, n_foc)), rbind(ref, foc))
  names(dat) <- c("x", paste0("i", seq_along(b_ref)))
  dat
}

## Read a 2011 data file (dif_responseiteration*.csv or DIF_data.dat):
## no header, first column = group, then 10 items; comma or tab separated.
read_dif_file <- function(file, items = ITEMS) {
  first <- readLines(file, n = 1)
  sep <- if (grepl(",", first)) "," else ""
  read.table(file, header = FALSE, sep = sep, col.names = c("x", items))
}

## ---- 2. MIMIC DIF (replaces Mplus diff_ab1.inp) -------------------------------

## lavaan translation of original/diff_ab1.inp:
##   f1 BY i1* i2-i10;  f1@1;  i1-i10 ON x;   (CATEGORICAL, WLSMV, delta)
mimic_syntax <- function(items = ITEMS, group = "x") {
  paste(c(
    paste0("f1 =~ NA*", paste(items, collapse = " + ")),   # f1 BY i1* i2-i10
    "f1 ~~ 1*f1",                                          # f1@1
    paste0(items, " ~ ", group)                            # i ON x (DIF effects)
  ), collapse = "\n")
}

fit_mimic <- function(data, items = ITEMS, group = "x") {
  lavaan::sem(mimic_syntax(items, group), data = data, ordered = items,
              estimator = "WLSMV", parameterization = "delta")
}

## The direct effects i ON x: est, se, z (= Mplus "Est./S.E."), p.
mimic_dif_effects <- function(fit, group = "x") {
  pe <- lavaan::parameterEstimates(fit)
  pe <- pe[pe$op == "~" & pe$rhs == group, ]
  data.frame(item = pe$lhs, est = pe$est, se = pe$se, z = pe$z, p = pe$pvalue,
             row.names = NULL)
}

## One MIMIC DIF analysis. The flag rule is the 2011 one,
##   Pseudo_Z <- (Est./S.E. >= 1.96),
## i.e. one-sided (alpha = .025), flagging only DIF favouring the focal group.
## A failed or non-converged fit returns NA instead of stopping the simulation.
run_mimic <- function(data, items = ITEMS, crit = 1.96) {
  na <- rep(NA_real_, length(items))
  warn <- FALSE
  res <- tryCatch(withCallingHandlers({
    fit <- fit_mimic(data, items)
    if (!lavaan::lavInspect(fit, "converged")) stop("lavaan did not converge")
    mimic_dif_effects(fit)
  }, warning = function(w) { warn <<- TRUE; invokeRestart("muffleWarning") }),
  error = function(e) NULL)
  if (is.null(res)) return(list(z = na, flag = na, status = "error"))
  list(z = res$z, flag = res$z >= crit, status = if (warn) "warning" else "ok")
}

## ---- 3. IRT likelihood-ratio DIF (difR, as in the 2011 code) ------------------

## abc <- difLRT(DIF_data[,2:11], group = DIF_data[,1], focal.name = 1, purify = TRUE)
## LRTstat[i,] <- (abc$LRT >= abc$thr)
run_lrt <- function(data, items = ITEMS, purify = TRUE, alpha = 0.05) {
  na <- rep(NA_real_, length(items))
  res <- tryCatch(
    suppressWarnings(difR::difLRT(data[, items], group = data$x, focal.name = 1,
                                  alpha = alpha, purify = purify)),
    error = function(e) NULL)
  if (is.null(res)) return(list(stat = na, flag = na, status = "error"))
  list(stat = res$LRT, flag = res$LRT >= res$thr, status = "ok")
}

## ---- 4. Rasch tree DIF (psychotree; added in this rebuild) ---------------------

## Two steps, following Strobl, Kopf & Zeileis (2015) and Kopf, Zeileis & Strobl (2015):
##  1. raschtree(): a score-based test of whether the Rasch item parameters are the same
##     in both groups (with one binary covariate: a global DIF test, chi-square with
##     k - 1 df). If it rejects at alpha, the tree splits on group.
##  2. If the tree split, anchortest(): select anchor items with the iterative forward
##     class and MTT selection, then Wald-test every item against that anchor. The first
##     anchor item gets no test (it fixes the scale) and is counted as not flagged.
## If the tree does not split, no item is flagged.
## stat is the signed Wald z (positive = easier for the focal group); NA if not tested.
run_tree <- function(data, items = ITEMS, alpha = 0.05) {
  ## raschtree() calls partykit::mob() unqualified, so psychotree must be attached
  if (!"package:psychotree" %in% search()) suppressPackageStartupMessages(library(psychotree))
  na <- rep(NA_real_, length(items))
  dd <- data.frame(x = factor(data$x))
  dd$resp <- as.matrix(data[, items])
  res <- tryCatch(suppressWarnings({
    tree <- psychotree::raschtree(resp ~ x, data = dd, alpha = alpha)
    split <- length(partykit::nodeids(tree, terminal = TRUE)) > 1
    z <- p <- setNames(rep(NA_real_, length(items)), items)
    if (split) {
      at <- psychotools::anchortest(resp ~ x, data = dd, class = "forward", select = "MTT")
      tested <- sub("^resp", "", names(at$final_tests$test$tstat))
      z[tested] <- at$final_tests$test$tstat
      p[tested] <- at$final_tests$test$pvalues
    }
    list(split = split, z = unname(z), p = unname(p))
  }), error = function(e) NULL)
  if (is.null(res)) return(list(stat = na, flag = na, split = NA, status = "error"))
  flag <- if (res$split) !is.na(res$p) & res$p < alpha else rep(FALSE, length(items))
  list(stat = res$z, flag = flag, split = res$split, status = "ok")
}

## ---- 5. Reading Mplus 4.2 output (for checking lavaan against Mplus) ----------

## Parse the MODEL RESULTS block of an Mplus 4.2 .out file into lavaan's
## lhs / op / rhs naming. Only BY, ON, Thresholds and Variances are needed here.
read_mplus42_results <- function(file) {
  x <- readLines(file)
  x <- x[(grep("^MODEL RESULTS", x) + 1):(grep("^R-SQUARE", x) - 1)]
  rows <- list()
  section <- NA
  for (line in trimws(x)) {
    f <- strsplit(line, "\\s+")[[1]]
    if (line == "" || f[1] == "Estimates") next
    num <- suppressWarnings(as.numeric(f[-1]))
    if (length(f) == 4 && !anyNA(num)) {
      rows[[length(rows) + 1]] <- data.frame(section = section, param = f[1],
                                             est = num[1], se = num[2], z = num[3])
    } else {
      section <- line
    }
  }
  out <- do.call(rbind, rows)
  s <- strsplit(out$section, "\\s+")
  head <- tolower(vapply(s, `[`, "", 1))
  kind <- vapply(s, function(v) v[length(v)], "")
  par <- tolower(out$param)
  out$lhs <- ifelse(kind %in% c("BY", "ON"), head, sub("\\$.*", "", par))
  out$op  <- c(BY = "=~", ON = "~", Thresholds = "|", Variances = "~~")[kind]
  out$rhs <- ifelse(kind == "Thresholds", paste0("t", sub(".*\\$", "", par)), par)
  out[, c("lhs", "op", "rhs", "est", "se", "z")]
}
