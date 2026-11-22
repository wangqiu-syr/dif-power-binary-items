## =============================================================================
## 02_simulation.R
## Monte Carlo replication of the 2011 study (original/Mplus_DIFR_R.r), plus a
## third method. Per replication, simulate 1PL data for 50 reference + 50 focal
## examinees, then test every item for DIF with
##   - MIMIC (lavaan, WLSMV; was Mplus diff_ab1.inp)       flag: z >= 1.96
##   - IRT likelihood-ratio test (difR::difLRT, purified)   flag: LRT >= thr
##   - Rasch tree (psychotree) + anchored Wald tests        flag: tree splits and p < .05
##
## Run from the repository root, either from a terminal:
##   Rscript R/02_simulation.R              # 5 replications, as in 2011
##   Rscript R/02_simulation.R 1000         # 1000 replications
##   Rscript R/02_simulation.R 1000 12345   # ... with another seed
## or from the R / RStudio console (working directory = repository root):
##   n_rep <- 1000; source("R/02_simulation.R")
##   n_rep <- 1000; seed <- 12345; source("R/02_simulation.R")
##   n_rep <- 1000; methods <- "TREE"; source("R/02_simulation.R")   # one method only
##
## Each data set depends only on the seed and the replication number, not on
## which methods are run: the random-number state is saved after the data are
## generated and restored after the methods have run. So runs of different
## methods with the same seed analyse the same data sets.
## =============================================================================

source(file.path("R", "functions.R"))
dir.create("results", showWarnings = FALSE)

## Command-line arguments win; otherwise use n_rep / seed / methods set in the
## console; otherwise the defaults.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) n_rep <- as.integer(args[1])
if (length(args) >= 2) seed  <- as.integer(args[2])
if (!exists("n_rep",   inherits = FALSE)) n_rep   <- 5L          # 2011: rep1 <- 5
if (!exists("seed",    inherits = FALSE)) seed    <- 20110219L
if (!exists("methods", inherits = FALSE)) methods <- c("MIMIC", "LRT", "TREE")

## Skip any method whose package is not installed.
needs <- c(MIMIC = "lavaan", LRT = "difR", TREE = "psychotree")
missing_pkg <- methods[!vapply(needs[methods], requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkg)) {
  warning("Skipping ", paste(missing_pkg, collapse = ", "), ": install.packages(c(",
          paste0('"', needs[missing_pkg], '"', collapse = ", "), "))", call. = FALSE)
  methods <- setdiff(methods, missing_pkg)
}
cat(sprintf("Running %d replications, seed %d, methods: %s\n",
            n_rep, seed, paste(methods, collapse = ", ")))

## ---- Design (2011 values) -------------------------------------------------------
n_ref <- 50
n_foc <- 50
b_ref <- c(-1,   -1, -1, 0,  0,  0, 0, 1, 1, 1)    # rdif
b_foc <- c(-1.8, -1, -1, 0, -.8, 0, 0, 1, 1, .2)   # fdif
dif_items <- which(b_ref != b_foc)                 # i1, i5, i10: focal easier by 0.8
k <- length(b_ref)

## ---- Replications --------------------------------------------------------------------
empty <- function() matrix(NA_real_, n_rep, k, dimnames = list(NULL, ITEMS))
stat   <- list(MIMIC = empty(), LRT = empty(), TREE = empty())  # z, LRT, z
flags  <- list(MIMIC = empty(), LRT = empty(), TREE = empty())  # 2011 names: Pseudo_Z, LRTstat
status <- list(MIMIC = rep("skipped", n_rep), LRT = rep("skipped", n_rep),
               TREE = rep("skipped", n_rep))
tree_split <- rep(NA, n_rep)

set.seed(seed)
t0 <- Sys.time()
for (r in seq_len(n_rep)) {
  dat <- simulate_dif_data(b_ref, b_foc, n_ref, n_foc)
  rng <- .Random.seed                       # data stream for the next replication

  if ("MIMIC" %in% methods) {
    m <- run_mimic(dat)
    stat$MIMIC[r, ] <- m$z; flags$MIMIC[r, ] <- m$flag; status$MIMIC[r] <- m$status
  }
  if ("LRT" %in% methods) {
    l <- run_lrt(dat)
    stat$LRT[r, ] <- l$stat; flags$LRT[r, ] <- l$flag; status$LRT[r] <- l$status
  }
  if ("TREE" %in% methods) {
    set.seed(seed + 100000L + r)            # anchortest() breaks ties at random
    tr <- run_tree(dat)
    stat$TREE[r, ] <- tr$stat; flags$TREE[r, ] <- tr$flag; status$TREE[r] <- tr$status
    tree_split[r] <- tr$split
  }

  assign(".Random.seed", rng, envir = globalenv())   # methods do not shift the data
  if (r %% max(1, n_rep %/% 10) == 0 || r == n_rep) {
    cat(sprintf("replication %d / %d  (%.0f s)\n", r, n_rep,
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
Pseudo_Z <- flags$MIMIC     # 2011 names
LRTstat  <- flags$LRT

## ---- Results ---------------------------------------------------------------------------
if ("MIMIC" %in% methods) { cat("\nPseudo_Z (MIMIC: 1 = flagged as DIF)\n"); print(head(Pseudo_Z, 10)) }
if ("LRT"   %in% methods) { cat("\nLRTstat (difLRT: 1 = flagged as DIF)\n"); print(head(LRTstat, 10)) }
if ("TREE"  %in% methods) { cat("\nRasch tree (1 = flagged as DIF)\n"); print(head(flags$TREE, 10)) }
if (n_rep > 10) cat("(first 10 replications shown)\n")

cat("\n")
for (mt in methods) {
  tb <- table(status[[mt]])
  cat(sprintf("%-5s fits: %s\n", mt, paste(names(tb), tb, sep = " = ", collapse = ", ")))
}
if ("TREE" %in% methods) {
  cat(sprintf("Rasch tree split on group (global DIF detected) in %.1f%% of replications\n",
              100 * mean(tree_split, na.rm = TRUE)))
}

## Flag rate per item = power for DIF items, Type I error rate for the others.
## Failed fits (NA) are left out of the denominator.
labels <- c(MIMIC = "MIMIC (z >= 1.96)", LRT = "LRT (difLRT, purified)",
            TREE = "Rasch tree + anchored Wald")
by_item <- data.frame(item = ITEMS, b_ref = b_ref, b_foc = b_foc,
                      dif = seq_len(k) %in% dif_items)
for (mt in methods) by_item[[mt]] <- round(colMeans(flags[[mt]], na.rm = TRUE), 3)
summary_tab <- data.frame(
  method = unname(labels[methods]),
  power = vapply(methods, function(mt) mean(flags[[mt]][, dif_items], na.rm = TRUE), 0),
  type1 = vapply(methods, function(mt) mean(flags[[mt]][, -dif_items], na.rm = TRUE), 0),
  row.names = NULL)
summary_tab[, 2:3] <- round(summary_tab[, 2:3], 3)

cat(sprintf("\nFlag rates over %d replications (seed %d)\n", n_rep, seed))
print(by_item, row.names = FALSE)
cat("\n"); print(summary_tab, row.names = FALSE)

## Long file with every statistic, so other flag rules can be applied later
## (e.g. two-sided |z| >= 1.96 for MIMIC).
long <- do.call(rbind, lapply(methods, function(mt) data.frame(
  rep = rep(seq_len(n_rep), k), item = rep(ITEMS, each = n_rep), method = mt,
  stat = c(stat[[mt]]), flag = c(flags[[mt]]))))
long <- long[order(long$rep, long$method, match(long$item, ITEMS)), ]
if ("TREE" %in% methods) {
  long$tree_split <- ifelse(long$method == "TREE", tree_split[long$rep], NA)
}

tag <- sprintf("rep%d_seed%d_%s", n_rep, seed, paste(methods, collapse = "-"))
write.csv(long, file.path("results", sprintf("sim_replications_%s.csv", tag)), row.names = FALSE)
write.csv(by_item, file.path("results", sprintf("sim_item_rates_%s.csv", tag)), row.names = FALSE)
write.csv(summary_tab, file.path("results", sprintf("sim_summary_%s.csv", tag)), row.names = FALSE)
cat("\nWrote results/sim_*_", tag, ".csv\n", sep = "")
