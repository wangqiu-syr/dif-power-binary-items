## =============================================================================
## 04_power_pilot.R
## Power pilot (docs/power.md, section 5): every method's variants on the same data sets, to find
## where power is lost. Designs:
##   study      the 2011 design (50 + 50; DIF of 0.8 on i1, i5, i10, all favouring the focal group).
##              Seeds 101 and 102, replications 1-100: the same data sets as the 1,000-replication run.
##   null       no DIF at all (b_foc = b_ref), 50 + 50; seed 301
##   study_n150 the 2011 DIF, 150 + 150; seed 401; no glmer LRT (slow), everything else
## One file per replication in results/power/runs/ (not tracked); reruns skip finished files.
##   Rscript R/04_power_pilot.R [cores] [reps] [designs]
##   Rscript R/04_power_pilot.R 4 200 study,null,study_n150
## =============================================================================
source(file.path("R", "functions.R"))
source(file.path("R", "power", "lrt_variants.R"))
source(file.path("R", "power", "pilot_methods.R"))
suppressMessages({ library(lavaan); library(psychotree) })
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 4L
reps <- if (length(args) >= 2) as.integer(args[2]) else 200L
designs <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else c("study", "null", "study_n150")

b_ref <- c(-1, -1, -1, 0, 0, 0, 0, 1, 1, 1)
b_dif <- c(-1.8, -1, -1, 0, -.8, 0, 0, 1, 1, .2)
DESIGNS <- list(
  study      = list(b_foc = b_dif, n = 50L,  seeds = c(101L, 102L), lrt = TRUE),
  null       = list(b_foc = b_ref, n = 50L,  seeds = 301L,          lrt = TRUE),
  study_n150 = list(b_foc = b_dif, n = 150L, seeds = 401L,          lrt = FALSE))
out_dir <- file.path("results", "power", "runs"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Data sets in the 02_simulation.R order: set.seed(seed), then replication r is the r-th draw.
jobs <- list()
for (dn in designs) {
  D <- DESIGNS[[dn]]; per_seed <- ceiling(reps / length(D$seeds))
  for (s in D$seeds) {
    set.seed(s)
    for (r in seq_len(per_seed)) {
      sim <- simulate_dif_theta(b_ref, D$b_foc, D$n, D$n)
      jobs[[length(jobs) + 1]] <- list(design = dn, seed = s, rep = r, sim = sim, D = D)
    }
  }
}

## Methods added after the first pilot pass (fast; computed for existing files too)
extra_methods <- function(resp, x, clean) {
  c(list(mh_none_nocc = mh_test(resp, x, correct = FALSE),
         mh_purified_nocc = mh_purified(resp, x, correct = FALSE),
         mh_oracle_nocc = mh_test(resp, x, function(j) setdiff(clean, j), correct = FALSE)),
    cml_variants(resp, x, clean))
}
to_long <- function(res, job, root_p, secs) {
  long <- do.call(rbind, lapply(names(res), function(m) data.frame(method = m, item = ITEMS,
                                                                  stat = res[[m]]$z, p = res[[m]]$p)))
  clean <- which(b_ref == job$D$b_foc)
  data.frame(design = job$design, seed = job$seed, rep = job$rep, n = job$D$n,
             dif = long$item %in% ITEMS[-clean], long, tree_root_p = root_p, secs = secs)
}

one <- function(job) {
  f <- file.path(out_dir, sprintf("%s_%d_%04d.rds", job$design, job$seed, job$rep))
  D <- job$D; dat <- job$sim$dat; resp <- as.matrix(dat[, ITEMS]); x <- dat$x
  clean <- which(b_ref == D$b_foc)
  if (file.exists(f)) {                                       # add the later methods if missing
    old <- readRDS(f)
    if (!"cml_oracle" %in% old$method) {
      saveRDS(rbind(old, to_long(extra_methods(resp, x, clean), job, old$tree_root_p[1], NA_real_)), f)
    }
    return(invisible(f))
  }
  t0 <- Sys.time()
  safe <- function(expr) tryCatch(expr, error = function(e) list(z = rep(NA_real_, 10), p = rep(NA_real_, 10)))
  res <- list(
    oracle_theta = safe(oracle_theta(resp, x, job$sim$theta)),
    marginal = safe(marginal_test(resp, x)),
    mh_none = safe(mh_test(resp, x)),
    mh_purified = safe(mh_purified(resp, x)),
    mh_oracle = safe(mh_test(resp, x, function(j) setdiff(clean, j))))
  res <- c(res, mimic_variants(dat, ITEMS, clean))
  tv <- tree_variants(resp, x, ITEMS, job$seed + 100000L + job$rep)
  res$anchortest <- tv$anchortest
  if (D$lrt) {
    l <- tryCatch(lrt_all(resp, x, clean), error = function(e) NULL)
    for (v in c("none", "difR", "standard", "oracle"))
      res[[paste0("lrt_", v)]] <- if (is.null(l)) list(z = rep(NA_real_, 10), p = rep(NA_real_, 10)) else
        list(z = l[[v]], p = pchisq(l[[v]], 1, lower.tail = FALSE))      # z holds the LRT statistic
  }
  res <- c(res, extra_methods(resp, x, clean))
  saveRDS(to_long(res, job, tv$tree_root_p, as.numeric(Sys.time() - t0, units = "secs")), f)
  invisible(f)
}

cat(sprintf("Power pilot: %d jobs (%s), %d cores\n", length(jobs), paste(designs, collapse = ", "), cores))
t0 <- Sys.time()
## big jobs first
ord <- order(-vapply(jobs, function(j) j$D$lrt, TRUE))
invisible(parallel::mclapply(jobs[ord], function(j) {
  set.seed(j$seed * 1000L + j$rep); tryCatch(one(j), error = function(e) message("job failed: ", conditionMessage(e)))
}, mc.cores = cores, mc.preschedule = FALSE))
cat(sprintf("Done in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
print(sessionInfo())
