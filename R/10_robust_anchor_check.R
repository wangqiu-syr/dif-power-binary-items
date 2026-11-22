## =============================================================================
## 10_robust_anchor_check.R
## Why the robust-anchor CML LRT exceeds .05 under one-way DIF of 0.5 logits at n = 250
## (docs/power_study.md, section 5). Two checks:
##  1) the study's own data sets (cell 25, replications 1-60, same seeds as R/08_power_study.R):
##     how often each item enters the anchor for clean item i2, and how far median(d) is pulled
##     toward the DIF items;
##  2) a normal approximation, d_i ~ N(shift_i, s^2) independently, anchor = the 5 (or 3) items
##     closest to the median: the Type I error that the selection alone produces.
##   Rscript R/10_robust_anchor_check.R     # writes results/power_study/robust_anchor_check.txt
## Seeds: as in R/08_power_study.R for part 1; set.seed(2) for part 2.
## =============================================================================
out_file <- file.path("results", "power_study", "robust_anchor_check.txt")
sink(out_file)
source(file.path("R", "08_power_study.R"))
cells <- study_cells(); cell <- cells[cells$cell == 25, ]
set.seed(1); nrep <- 60; inA <- matrix(0, nrep, 10); medpos <- numeric(nrep)
for (r in seq_len(nrep)) {
  seed <- MASTER_SEED + 1000L * 25 + r; set.seed(seed); sim <- simulate_cell(cell)
  resp <- as.matrix(sim$dat[, ITEMS]); x <- sim$dat$x
  d <- robust_anchor_d(resp, x); dist <- abs(d - median(d))
  ## anchor for clean item 2
  a <- setdiff(order(dist), 2)[1:5]; inA[r, a] <- 1
  medpos[r] <- median(d) - mean(d[-c(1, 5, 10)])
}
cat("share of replications in which each item is in the anchor for item 2 (n = 250, one-way 0.5):\n")
print(round(colMeans(inA), 2))
cat("median(d) minus the mean of the clean items' d (positive = pulled toward the DIF items):",
    round(mean(medpos), 3), "\n")
## 2) Normal approximation: d_i ~ N(shift_i, s^2); anchor = 5 nearest the median (excluding j)
approx_t1 <- function(dif, s, k = 5, R = 2e5) {
  sh <- c(dif, 0, 0, 0, dif, 0, 0, 0, 0, dif); D <- matrix(rnorm(R * 10, rep(sh, each = R), s), R)
  med <- apply(D, 1, median); out <- numeric(0)
  for (j in c(2, 3, 4, 6, 7, 8, 9)) {
    dist <- abs(D - med); dist[, j] <- Inf
    A <- t(apply(dist, 1, function(v) order(v)[1:k]))
    mA <- rowMeans(matrix(D[cbind(rep(seq_len(R), k), as.vector(A))], R))
    est <- D[, j] - mA; out <- c(out, mean(abs(est) / (s * sqrt(1 + 1 / k)) > qnorm(.975)))
  }
  mean(out)
}
set.seed(2)
for (s in c(.35, .25, .2)) cat(sprintf("normal approx, s = %.2f: Type I k=5 %.3f (dif .5), %.3f (dif .8); k=3 %.3f (dif .5)\n",
  s, approx_t1(.5, s, R = 4e4), approx_t1(.8, s, R = 4e4), approx_t1(.5, s, k = 3, R = 4e4)))
sink()
cat("Wrote", out_file, "\n")
