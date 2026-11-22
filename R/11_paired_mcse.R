# MCSE of the difference between two methods' rates in one Simulation 2 condition.
# The methods are applied to the same data sets, so the difference is paired by replication.
# For each method, replication r contributes (f_r - R t_r) / tbar (the ratio-estimator
# linearization of Appendix A2); the paired MCSE is the SD of the difference of these
# contributions over replications, divided by sqrt(n_rep).
# Output: results/power_study/paired_diff_mcse.csv (one row per condition, method pair, measure)
#         results/power_study/paired_diff_mcse_summary.txt
# Run from the repository root: Rscript R/11_paired_mcse.R

x <- readRDS("results/power_study/replications.rds")
rates <- read.csv("results/power_study/rates.csv")
cells <- unique(rates[, c("cell", "n", "impact", "direction", "dif")])

methods <- c("oracle_theta", "cml_oracle", "cml_purified", "mh_purified", "cml_robust5",
             "mimic2011", "mimic_msp", "tree_gated", "anchortest")
x <- x[x$method %in% methods & x$item != "bundle" & !is.na(x$p), ]
x$flag <- x$p < .05
dif_items <- c("i1", "i5", "i10")

# Per replication, method and measure: flags f and tests t
contrib <- function(d) {
  f <- tapply(d$flag, d$rep, sum)
  t <- tapply(d$flag, d$rep, length)
  R <- sum(f) / sum(t)
  u <- (f - R * t) / mean(t)
  data.frame(rep = as.integer(names(u)), u = as.numeric(u), R = R)
}

out <- list()
for (i in seq_len(nrow(cells))) {
  cc <- cells[i, ]
  xc <- x[x$cell == cc$cell, ]
  dif_cell <- cc$direction != "none"
  for (measure in c("power", "type1")) {
    if (measure == "power" && !dif_cell) next
    sel <- if (measure == "power") xc$item %in% dif_items else
      if (dif_cell) !(xc$item %in% dif_items) else rep(TRUE, nrow(xc))
    xm <- xc[sel, ]
    cl <- lapply(split(xm, droplevels(xm$method)), contrib)
    ms <- intersect(methods, names(cl))
    for (a in seq_along(ms)) for (b in seq_along(ms)) if (a < b) {
      m <- merge(cl[[ms[a]]], cl[[ms[b]]], by = "rep")
      dd <- m$u.x - m$u.y
      out[[length(out) + 1]] <- data.frame(
        cell = cc$cell, n = cc$n, impact = cc$impact, direction = cc$direction, dif = cc$dif,
        measure = measure, method_a = ms[a], method_b = ms[b],
        rate_a = m$R.x[1], rate_b = m$R.y[1], diff = m$R.x[1] - m$R.y[1],
        mcse_paired = sd(dd) / sqrt(nrow(m)),
        mcse_unpaired = sqrt(var(m$u.x) + var(m$u.y)) / sqrt(nrow(m)))
    }
  }
}
res <- do.call(rbind, out)
write.csv(res, "results/power_study/paired_diff_mcse.csv", row.names = FALSE)

# Summary: among the realistic methods (no benchmarks) and against the clean-anchor benchmark
real <- c("cml_purified", "mh_purified", "cml_robust5", "mimic2011", "mimic_msp",
          "tree_gated", "anchortest")
sink("results/power_study/paired_diff_mcse_summary.txt")
cat("Paired MCSE of the difference between two methods' rates, Simulation 2\n\n")
for (measure in c("power", "type1")) {
  r <- res[res$measure == measure, ]
  rr <- r[r$method_a %in% c(real, "cml_oracle") & r$method_b %in% c(real, "cml_oracle"), ]
  pur <- c("cml_purified", "mh_purified", "mimic_msp")
  rp <- r[r$method_a %in% pur & r$method_b %in% pur, ]
  cat(sprintf("%s, all pairs of realistic methods and the clean anchor: median %.3f, max %.3f; unpaired max %.3f\n",
              measure, median(rr$mcse_paired), max(rr$mcse_paired), max(rr$mcse_unpaired)))
  cat(sprintf("%s, pairs among the three purified methods: median %.3f, max %.3f\n",
              measure, median(rp$mcse_paired), max(rp$mcse_paired)))
}
sink()
cat(readLines("results/power_study/paired_diff_mcse_summary.txt"), sep = "\n")
