## =============================================================================
## 07_lrt_calibration.R
## Is the GLMM likelihood-ratio test itself calibrated at 50 + 50 with 10 binary items?
## (docs/power.md, section 5.3). Data with NO DIF and no impact; every item tested once against all
## other items (no purification, so no contamination and no accumulation). Four ways to compute
## the same Rasch LRT:
##   difR_deviance  difR's statistic: difference of lme4's deviance(), Laplace (bobyqa); for the first
##                  20 replications also with lme4's default optimizer, as difR calls it
##   laplace_lr     the likelihood ratio from the same fits: 2 x difference of logLik()
##   agq15_lr       the likelihood ratio with adaptive Gauss-Hermite quadrature (15 nodes)
##   cml_lr         conditional ML (psychotools::raschmodel): item j split into a reference and a
##                  focal "virtual item"; no assumption about the theta distribution
##   Rscript R/07_lrt_calibration.R [cores] [reps]
## =============================================================================
source(file.path("R", "functions.R")); source(file.path("R", "power", "lrt_variants.R"))
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 4L
reps <- if (length(args) >= 2) as.integer(args[2]) else 60L
b <- c(-1, -1, -1, 0, 0, 0, 0, 1, 1, 1); seed <- 501L
out_dir <- file.path("results", "power", "runs"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## One set of fits gives both statistics: difR's (difference of lme4 deviance(), which for a GLMM is
## the conditional deviance, not -2 log-likelihood) and the likelihood ratio (-2 logLik difference).
glmm_lrt <- function(resp, group, ctl, nAGQ = 1) {
  d <- lrt_long(resp, group)
  m0 <- lme4::glmer(y ~ item - 1 + member + (1 | person), data = d, family = binomial, control = ctl, nAGQ = nAGQ)
  t(vapply(seq_len(ncol(resp)), function(j) {
    d$X <- as.numeric(d$item == j)
    m1 <- lme4::glmer(y ~ item - 1 + member + member:X + (1 | person), data = d, family = binomial,
                      control = ctl, nAGQ = nAGQ)
    c(dev = stats::deviance(m0) - stats::deviance(m1),
      lr = 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))) }, c(dev = 0, lr = 0)))
}
cml_lrt <- function(resp, group) {
  ll0 <- as.numeric(logLik(psychotools::raschmodel(resp)))
  vapply(seq_len(ncol(resp)), function(j) {
    r2 <- cbind(resp[, -j], jR = ifelse(group == 0, resp[, j], NA), jF = ifelse(group == 1, resp[, j], NA))
    2 * (as.numeric(logLik(psychotools::raschmodel(r2))) - ll0) }, 0)
}

set.seed(seed); sets <- lapply(seq_len(reps), function(r) simulate_dif_data(b, b))
res <- parallel::mclapply(seq_len(reps), function(r) {
  f <- file.path(out_dir, sprintf("calib_%d_%04d.rds", seed, r))
  if (file.exists(f)) return(readRDS(f))
  dat <- sets[[r]]; resp <- as.matrix(dat[, ITEMS]); g <- dat$x
  fast <- lme4::glmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
  la <- glmm_lrt(resp, g, fast); aq <- glmm_lrt(resp, g, fast, nAGQ = 15)
  s <- list(difR_deviance = la[, "dev"], laplace_lr = la[, "lr"], agq15_lr = aq[, "lr"], cml_lr = cml_lrt(resp, g),
            difR_deviance_default_optimizer = if (r <= 20) glmm_lrt(resp, g, lme4::glmerControl())[, "dev"] else NA_real_)
  out <- data.frame(rep = r, item = ITEMS, s)
  saveRDS(out, f); out
}, mc.cores = cores, mc.preschedule = FALSE)
x <- do.call(rbind, res)
eng <- c("difR_deviance", "difR_deviance_default_optimizer", "laplace_lr", "agq15_lr", "cml_lr")
tab <- data.frame(engine = eng,
                  tests = vapply(eng, function(e) sum(!is.na(x[[e]])), 0),
                  mean_stat = vapply(eng, function(e) mean(x[[e]], na.rm = TRUE), 0),
                  type1 = vapply(eng, function(e) mean(x[[e]] > qchisq(.95, 1), na.rm = TRUE), 0),
                  q95 = vapply(eng, function(e) unname(quantile(x[[e]], .95, na.rm = TRUE)), 0))
tab$type1_se <- sqrt(tab$type1 * (1 - tab$type1) / tab$tests)
write.csv(transform(tab, mean_stat = round(mean_stat, 3), type1 = round(type1, 3), q95 = round(q95, 2),
                    type1_se = round(type1_se, 3)),
          file.path("results", "power", "lrt_calibration.csv"), row.names = FALSE)
print(tab, row.names = FALSE)
cat("chi-square(1): mean 1, 95th percentile 3.84\n")
print(sessionInfo())
