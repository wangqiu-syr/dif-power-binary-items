## Simulation 1: regenerate replications 1-2 of seed 101 and compare with the stored statistics.
## The stored run used R 4.5.1 and lavaan 0.7.2 on Windows; small differences in the last digits
## of the MIMIC z are expected with other versions, so the check allows 0.01.
if (!skip_if_missing(c("lavaan", "psychotree", "psychotools", "partykit"))) {
  source(file.path("R", "functions.R"))
  stored <- read.csv(file.path("results", "sim_replications_rep1000_MIMIC-LRT-TREE.csv"))
  b_ref <- c(-1, -1, -1, 0, 0, 0, 0, 1, 1, 1); b_foc <- c(-1.8, -1, -1, 0, -0.8, 0, 0, 1, 1, 0.2)
  seed <- 101L; set.seed(seed)
  has_difR <- requireNamespace("difR", quietly = TRUE)
  for (r in 1:2) {
    dat <- simulate_dif_data(b_ref, b_foc); rng <- .Random.seed
    st <- function(m) stored[stored$seed == seed & stored$rep == r & stored$method == m, ][match(ITEMS, stored$item[stored$seed == seed & stored$rep == r & stored$method == m]), ]
    m <- run_mimic(dat)
    check(sprintf("replication %d: MIMIC z within 0.01 of the stored run", r), max(abs(m$z - st("MIMIC")$stat)) < 0.01)
    check(sprintf("replication %d: MIMIC flags as stored", r), all(as.integer(m$flag) == st("MIMIC")$flag))
    set.seed(seed + 100000L + r); tr <- run_tree(dat)
    s <- st("TREE")
    check(sprintf("replication %d: Rasch tree split and flags as stored", r),
          tr$split == as.logical(s$tree_split[1]) && all(as.integer(tr$flag) == s$flag))
    check(sprintf("replication %d: Rasch tree z within 0.01 of the stored run", r),
          isTRUE(all.equal(is.na(tr$stat), is.na(s$stat))) && max(c(0, abs(tr$stat - s$stat)), na.rm = TRUE) < 0.01)
    if (has_difR && r == 1) {
      l <- run_lrt(dat)
      check("replication 1: difLRT statistics within 0.05 of the stored run", max(abs(l$stat - st("LRT")$stat)) < 0.05)
    }
    assign(".Random.seed", rng, envir = globalenv())
  }
  if (!has_difR) { .results$n["skip"] <- .results$n["skip"] + 1L; cat("  [SKIP] difLRT check needs difR\n") }
}
