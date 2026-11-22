## =============================================================================
## 08_power_study.R
## The full power study (docs/power_study.md; planned in docs/power.md, section 7).
## 40 cells: examinees per group (50, 100, 150, 250) x DIF direction (one-way, mixed) x impact
## (0, focal group 0.5 SD lower) x DIF size (0.5, 0.8 logits), plus 8 cells with no DIF
## (n x impact). DIF is on i1, i5 and i10; "one-way" makes all three easier for the focal group,
## "mixed" makes i5 harder for it. Rasch items, b = (-1, -1, -1, 0, 0, 0, 0, 1, 1, 1).
## One file per replication in results/power/study_runs/ (not tracked); reruns skip finished files.
##   Rscript R/08_power_study.R [cores] [reps] [cells] [out_dir]
##   Rscript R/08_power_study.R 4 300                 # all cells
##   Rscript R/08_power_study.R 4 20 all results/power/study_pilot_runs   # the pilot
## Seeds: replication r of cell c uses set.seed(20261119 + 1000 * c + r).
## =============================================================================
source(file.path("R", "functions.R"))
source(file.path("R", "power", "lrt_variants.R"))
source(file.path("R", "power", "pilot_methods.R"))
source(file.path("R", "power", "study_methods.R"))
suppressMessages({ library(lavaan); library(psychotree) })
args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1) as.integer(args[1]) else 4L
reps <- if (length(args) >= 2) as.integer(args[2]) else 300L
which_cells <- if (length(args) >= 3 && args[3] != "all") as.integer(strsplit(args[3], ",")[[1]]) else NULL
out_dir <- if (length(args) >= 4) args[4] else file.path("results", "power", "study_runs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

B_REF <- c(-1, -1, -1, 0, 0, 0, 0, 1, 1, 1)
DIF_ITEMS <- c(1L, 5L, 10L)
MASTER_SEED <- 20261119L

study_cells <- function() {
  g <- expand.grid(dif = c(0.5, 0.8), impact = c(0, 0.5), direction = c("one-way", "mixed"),
                   n = c(50L, 100L, 150L, 250L), stringsAsFactors = FALSE)
  z <- expand.grid(dif = 0, impact = c(0, 0.5), direction = "none", n = c(50L, 100L, 150L, 250L),
                   stringsAsFactors = FALSE)
  cells <- rbind(g, z)[, c("n", "direction", "impact", "dif")]
  cells$cell <- seq_len(nrow(cells))
  ## difLRT (glmer, slow) only where it is the 2011 reference: the 2011 design and its no-DIF twin
  cells$diflrt <- cells$n == 50 & cells$impact == 0 &
    ((cells$direction == "one-way" & cells$dif == 0.8) | cells$direction == "none")
  cells
}
true_shift <- function(cell) {             # + = easier for the focal group, per item
  s <- numeric(length(B_REF))
  if (cell$direction == "one-way") s[DIF_ITEMS] <- cell$dif
  if (cell$direction == "mixed") s[DIF_ITEMS] <- cell$dif * c(1, -1, 1)
  s
}
simulate_cell <- function(cell) {
  b_foc <- B_REF - true_shift(cell)
  th_r <- rnorm(cell$n); ref <- respon_1pl(B_REF, th_r)
  th_f <- rnorm(cell$n, -cell$impact); foc <- respon_1pl(b_foc, th_f)
  dat <- data.frame(c(rep(0L, cell$n), rep(1L, cell$n)), rbind(ref, foc))
  names(dat) <- c("x", ITEMS)
  list(dat = dat, theta = c(th_r, th_f))
}

one_rep <- function(cell, r, out_dir) {
  f <- file.path(out_dir, sprintf("cell%02d_%04d.rds", cell$cell, r))
  if (file.exists(f)) return(invisible(f))
  seed <- MASTER_SEED + 1000L * cell$cell + r
  set.seed(seed); sim <- simulate_cell(cell)
  t0 <- Sys.time(); dat <- sim$dat; resp <- as.matrix(dat[, ITEMS]); x <- dat$x
  shift <- true_shift(cell); clean <- which(shift == 0)
  safe <- function(expr) tryCatch(expr, error = function(e) signed(rep(NA_real_, 10)))
  res <- list(
    oracle_theta = safe(oracle_theta(resp, x, sim$theta)),
    mimic2011 = safe(mimic_2011(dat, ITEMS)),
    mh_purified = safe(mh_purified(resp, x, correct = FALSE)),
    cml_purified = safe(cml_purified_test(resp, x)),
    cml_robust5 = safe(cml_robust_test(resp, x, k = 5)),
    cml_robust3 = safe(cml_robust_test(resp, x, k = 3)),
    cml_oracle = safe(cml_test(resp, x, function(j) setdiff(if (length(clean) == 10) seq_len(10) else clean, j))))
  msp <- tryCatch(mimic_msp(dat, ITEMS), error = function(e) c(signed(rep(NA_real_, 10)), passes = NA))
  res$mimic_msp <- msp
  set.seed(seed + 7L)                                         # anchortest() breaks ties at random
  tv <- tryCatch(tree_variants(resp, x, ITEMS, seed + 7L), error = function(e) NULL)
  res$anchortest <- if (is.null(tv)) signed(rep(NA_real_, 10)) else signed(tv$anchortest$z)
  if (!is.null(tv)) res$anchortest$p <- tv$anchortest$p
  if (cell$diflrt) {
    s <- tryCatch(lrt_variant(resp, x, "difR"), error = function(e) rep(NA_real_, 10))
    res$difLRT <- list(z = s, p = pchisq(s, 1, lower.tail = FALSE))      # unsigned: z holds the LRT
  }
  long <- do.call(rbind, lapply(names(res), function(m)
    data.frame(method = m, item = ITEMS, stat = res[[m]]$z, p = res[[m]]$p)))
  bu <- tryCatch(bundle_test(resp, x, DIF_ITEMS), error = function(e) list(z = NA_real_, p = NA_real_))
  long <- rbind(long, data.frame(method = "bundle", item = "bundle", stat = bu$z, p = bu$p))
  out <- data.frame(cell = cell$cell, rep = r, seed = seed, long,
                    true_shift = c(rep(shift, length(res)), sum(shift[DIF_ITEMS])),
                    tree_root_p = if (is.null(tv)) NA_real_ else tv$tree_root_p,
                    msp_passes = msp$passes, robust5_fallback = isTRUE(res$cml_robust5$fallback),
                    secs = as.numeric(Sys.time() - t0, units = "secs"))
  saveRDS(out, f); invisible(f)
}

if (sys.nframe() == 0L) {
  cells <- study_cells()
  run <- if (is.null(which_cells)) cells else cells[cells$cell %in% which_cells, ]
  jobs <- do.call(c, lapply(seq_len(nrow(run)), function(i) lapply(seq_len(reps), function(r) list(cell = run[i, ], r = r))))
  ## slow jobs first: difLRT cells, then by n
  key <- vapply(jobs, function(j) j$cell$diflrt * 1e4 + j$cell$n, 0)
  jobs <- jobs[order(-key)]
  cat(sprintf("Power study: %d cells x %d replications = %d jobs, %d cores -> %s\n",
              nrow(run), reps, length(jobs), cores, out_dir))
  t0 <- Sys.time()
  invisible(parallel::mclapply(jobs, function(j) tryCatch(one_rep(j$cell, j$r, out_dir),
                                                          error = function(e) message("job failed: ", conditionMessage(e))),
                               mc.cores = cores, mc.preschedule = FALSE))
  cat(sprintf("Done in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
  print(sessionInfo())
}
