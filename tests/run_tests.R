## =============================================================================
## run_tests.R
## Runs every tests/test_*.R from the repository root and stops with an error if any check fails.
##   Rscript tests/run_tests.R
## Tests that need a package that is not installed are reported as skipped, not failed.
## Base R only (no testing framework needed).
## =============================================================================
if (!file.exists(file.path("R", "functions.R"))) stop("Run from the repository root: Rscript tests/run_tests.R")
if (!isTRUE(l10n_info()[["UTF-8"]])) invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))
.results <- new.env()
.results$n <- c(pass = 0L, fail = 0L, skip = 0L)
check <- function(desc, ok) {
  ok <- isTRUE(all(ok))
  .results$n[if (ok) "pass" else "fail"] <- .results$n[if (ok) "pass" else "fail"] + 1L
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", desc))
}
skip_if_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    .results$n["skip"] <- .results$n["skip"] + 1L
    cat(sprintf("  [SKIP] needs %s\n", paste(miss, collapse = ", ")))
    return(TRUE)
  }
  FALSE
}
for (f in sort(list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE))) {
  cat(basename(f), "\n")
  t0 <- Sys.time()
  tryCatch(sys.source(f, envir = new.env(parent = globalenv())),
           error = function(e) check(paste("ran without error:", conditionMessage(e)), FALSE))
  cat(sprintf("  (%.0f s)\n", as.numeric(Sys.time() - t0, units = "secs")))
}
n <- .results$n
cat(sprintf("\n%d passed, %d failed, %d skipped\n", n[["pass"]], n[["fail"]], n[["skip"]]))
if (n[["fail"]] > 0) quit(status = 1)
