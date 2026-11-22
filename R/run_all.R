## =============================================================================
## run_all.R
## One entry point for everything in the report. Run from the repository root:
##   Rscript R/run_all.R                 # QUICK (default): a few minutes
##   Rscript R/run_all.R FULL [cores]    # FULL: every simulation again, about 12 hours on 8 cores
##   Rscript R/run_all.R FULL 8 sim1,study   # FULL, only some steps (sim1, pilot, calibration, study)
##
## QUICK rebuilds the report's tables (R/report_tables.R) and figures (R/report_figures.R) from the
## stored results, and runs the tests (tests/run_tests.R), which also regenerate one replication of
## each simulation and compare it with the stored statistics. Nothing in results/ is overwritten
## except the three files that R/report_tables.R writes, which it reproduces exactly.
##
## FULL also reruns the simulations, in this order, and overwrites their files in results/:
##   sim1         Simulation 1: R/02_simulation.R as 8 jobs of 125 replications, seeds 101-108,
##                in parallel; then R/03_figures.R merges them into the rep1000 files and draws
##                results/figures/fig1, fig2. The per-seed files are deleted after the merge.
##   pilot        the power pilot: R/04_power_pilot.R (2011 design and no-DIF design, 100
##                replications each; 150 per group, 200 replications), then R/06_power_summary.R
##   calibration  the LRT calibration: R/07_lrt_calibration.R, 60 replications
##   study        Simulation 2: R/08_power_study.R (40 cells x 300 replications), then
##                R/09_power_study_summary.R and R/10_robust_anchor_check.R
## and then everything QUICK does. Per-replication files go to results/power/runs/ and
## results/power/study_runs/ (not tracked); a rerun skips replications whose file exists, so an
## interrupted FULL run can simply be restarted. Delete those folders to start from scratch.
##
## Not rerun: the exact power calculations (results/power/exact_*.csv). Their code is released
## with the follow-up study (README, "What is not in this repository").
## Seeds: every script sets its own (README, "Controlled randomness").
## =============================================================================

if (!file.exists(file.path("R", "functions.R"))) stop("Run from the repository root: Rscript R/run_all.R")
args  <- commandArgs(trailingOnly = TRUE)
mode  <- toupper(if (length(args) >= 1) args[1] else "QUICK")
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, min(8L, parallel::detectCores()))
steps <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else c("sim1", "pilot", "calibration", "study")
stopifnot(mode %in% c("QUICK", "FULL"), steps %in% c("sim1", "pilot", "calibration", "study"))

rscript <- file.path(R.home("bin"), "Rscript")
dir.create("logs", showWarnings = FALSE)
run <- function(script, ..., log = NULL) {
  a <- c(script, as.character(c(...)))
  cat(sprintf("[%s] Rscript %s\n", format(Sys.time(), "%H:%M:%S"), paste(a, collapse = " ")))
  t0 <- Sys.time()
  out <- if (is.null(log)) "" else file.path("logs", log)
  status <- system2(rscript, a, stdout = out, stderr = out)
  if (!identical(status, 0L)) stop(script, " failed (exit status ", status, ")", if (!is.null(log)) paste0("; see logs/", log))
  cat(sprintf("           done in %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
}

if (mode == "FULL") {
  need <- c("lavaan", "difR", "psychotree", "psychotools", "partykit", "lme4")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("FULL needs: install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))")
  cat(sprintf("FULL run on %d cores: %s\n", cores, paste(steps, collapse = ", ")))

  if ("sim1" %in% steps) {
    seeds <- 101:108
    old <- list.files("results", pattern = "^sim_(replications|item_rates|summary)_rep[0-9]+_seed[0-9]+_", full.names = TRUE)
    file.remove(old)                                          # stale per-seed files would be merged too
    job <- function(s) run(file.path("R", "02_simulation.R"), 125, s, log = sprintf("sim1_seed%d.log", s))
    if (.Platform$OS.type == "unix" && cores > 1) {
      res <- parallel::mclapply(seeds, function(s) tryCatch(job(s), error = function(e) conditionMessage(e)),
                                mc.cores = min(cores, length(seeds)), mc.preschedule = FALSE)
      bad <- vapply(res, is.character, TRUE)
      if (any(bad)) stop(paste(unlist(res[bad]), collapse = "\n"))
    } else {
      for (s in seeds) job(s)                                 # Windows: one job after another
    }
    run(file.path("R", "03_figures.R"), log = "sim1_merge.log")
    file.remove(list.files("results", pattern = "^sim_(replications|item_rates|summary)_rep125_seed[0-9]+_", full.names = TRUE))
  }
  if ("pilot" %in% steps) {
    run(file.path("R", "04_power_pilot.R"), cores, 100, "study,null", log = "pilot_n50.log")
    run(file.path("R", "04_power_pilot.R"), cores, 200, "study_n150", log = "pilot_n150.log")
    run(file.path("R", "06_power_summary.R"), log = "pilot_summary.log")
  }
  if ("calibration" %in% steps) run(file.path("R", "07_lrt_calibration.R"), cores, 60, log = "lrt_calibration.log")
  if ("study" %in% steps) {
    run(file.path("R", "08_power_study.R"), cores, 300, log = "study.log")
    run(file.path("R", "09_power_study_summary.R"), log = "study_summary.log")
    run(file.path("R", "10_robust_anchor_check.R"), log = "robust_anchor_check.log")
  }
}

## QUICK (and the end of FULL): tables, figures, tests
run(file.path("R", "report_tables.R"), log = "report_tables.log")
run(file.path("R", "report_figures.R"), log = "report_figures.log")
run(file.path("R", "11_paired_mcse.R"), log = "paired_mcse.log")
run(file.path("tests", "run_tests.R"))
cat("All done. Build the PDF with: bash report/pdf/build.sh\n")
