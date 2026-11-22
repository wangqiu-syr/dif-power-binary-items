## =============================================================================
## 06_power_summary.R
## Summaries of the power pilot (R/04_power_pilot.R; docs/power.md, section 5).
##   Rscript R/06_power_summary.R
## Reads results/power/runs/ (not tracked); writes results/power/pilot_rates.csv,
## pilot_item_rates.csv, pilot_replications.csv.gz and results/figures/fig3_power_pilot.png.
## =============================================================================
f <- list.files(file.path("results", "power", "runs"), pattern = "^(study|null).*\\.rds$", full.names = TRUE)
x <- do.call(rbind, lapply(f, readRDS))
out <- file.path("results", "power")
write.csv(x[, setdiff(names(x), "secs")], gzfile(file.path(out, "pilot_replications.csv.gz")), row.names = FALSE)
a <- .05
## The tree as in the 2011 rebuild: flagged only if the global test splits (root p < .05)
tg <- x[x$method == "anchortest", ]; tg$method <- "tree_gated"
tg$p <- ifelse(tg$tree_root_p < a, tg$p, 1)
x <- rbind(x, tg)
x$flag <- x$p < a
LABELS <- c(oracle_theta = "Benchmark: theta known", marginal = "Benchmark: proportions (no impact)",
  mimic2011 = "MIMIC as in 2011 (no f ~ x)", mimic_eqload = "MIMIC 2011, equal loadings",
  mimic_anchored = "MIMIC, f ~ x, all other items anchor", mimic_oracle = "MIMIC, f ~ x, true clean anchor",
  lrt_none = "LRT, no purification", lrt_difR = "LRT, difR purification (2011)",
  lrt_standard = "LRT, standard purification", lrt_oracle = "LRT, true clean anchor",
  mh_none = "Mantel-Haenszel, total score", mh_purified = "Mantel-Haenszel, purified",
  mh_oracle = "Mantel-Haenszel, true clean anchor",
  mh_none_nocc = "MH, no continuity correction", mh_purified_nocc = "MH, no correction, purified",
  mh_oracle_nocc = "MH, no correction, true clean anchor",
  cml_none = "CML LRT, all other items anchor", cml_purified = "CML LRT, standard purification",
  cml_oracle = "CML LRT, true clean anchor",
  tree_gated = "Rasch tree (gate) + anchortest (2011 rebuild)", anchortest = "anchortest, no gate")
rates <- do.call(rbind, lapply(split(x, list(x$design, x$method), drop = TRUE), function(z) {
  d <- z[z$dif, ]; c0 <- z[!z$dif, ]
  data.frame(design = z$design[1], n_per_group = z$n[1], method = z$method[1], label = LABELS[[z$method[1]]],
             reps = length(unique(paste(z$seed, z$rep))),
             power = if (nrow(d)) mean(d$flag, na.rm = TRUE) else NA_real_,
             type1 = mean(c0$flag, na.rm = TRUE), na = sum(is.na(z$p)))
}))
## size-adjusted power in the study design: critical p from the null design's clean items
null_q <- tapply(x$p[x$design == "null"], x$method[x$design == "null"], quantile, probs = a, na.rm = TRUE)
rates$power_size_adjusted <- NA_real_
for (i in which(rates$design == "study")) {
  m <- rates$method[i]; if (is.na(null_q[m]) || m == "tree_gated") next   # gated p is 1 when the tree does not split
  d <- x[x$design == "study" & x$method == m & x$dif, ]
  rates$power_size_adjusted[i] <- mean(d$p < null_q[m], na.rm = TRUE)
}
split_rate <- tapply(x$tree_root_p[x$method == "anchortest" & x$item == "i1"] < a,
                     x$design[x$method == "anchortest" & x$item == "i1"], mean)
ord <- names(LABELS)
rates <- rates[order(match(rates$design, c("study", "study_n150", "null")), match(rates$method, ord)), ]
item_rates <- aggregate(flag ~ design + method + item, x, mean)
rnd <- function(d) { i <- vapply(d, is.numeric, TRUE) & !names(d) %in% c("n_per_group", "reps", "na"); d[i] <- lapply(d[i], round, 3); d }
write.csv(rnd(rates), file.path(out, "pilot_rates.csv"), row.names = FALSE)
write.csv(rnd(item_rates), file.path(out, "pilot_item_rates.csv"), row.names = FALSE)
options(width = 200); print(rnd(rates)[, c("design", "method", "reps", "power", "power_size_adjusted", "type1", "na")], row.names = FALSE)
cat("\nRasch tree split rate (global test p < .05):\n"); print(round(split_rate, 3))

## Figure: power and Type I error, 2011 DIF, n = 50 and 150 per group
png(file.path("results", "figures", "fig3_power_pilot.png"), width = 2000, height = 1250, res = 170, type = "cairo")
show <- c("oracle_theta", "marginal", "mimic2011", "mimic_anchored", "mimic_oracle", "lrt_difR", "lrt_oracle",
          "cml_none", "cml_purified", "cml_oracle", "mh_purified_nocc", "mh_oracle_nocc", "tree_gated", "anchortest")
P <- sapply(c("study", "study_n150"), function(d) rates$power[match(paste(d, show), paste(rates$design, rates$method))])
T1 <- sapply(c("study", "study_n150"), function(d) rates$type1[match(paste(d, show), paste(rates$design, rates$method))])
cols <- c("#86b6ef", "#184f95"); lab <- unname(LABELS[show])
layout(matrix(1:2, 1), widths = c(1.55, 1))
par(mar = c(4.5, 17, 3, 1), las = 1, col.axis = "#3d3d3a", fg = "#8a8a85")
barplot(t(P[rev(seq_along(show)), ]), beside = TRUE, horiz = TRUE, xlim = c(0, 1), col = cols, border = NA,
        names.arg = rev(lab), cex.names = .78, main = "Power (i1, i5, i10)", col.main = "#1a1a19",
        xlab = "share of replications flagged")
abline(v = .8, lty = 3, col = "#8a8a85")
par(mar = c(4.5, 1, 3, 1))
barplot(t(T1[rev(seq_along(show)), ]), beside = TRUE, horiz = TRUE, xlim = c(0, max(.3, max(T1, na.rm = TRUE))),
        col = cols, border = NA, names.arg = rep("", length(show)), main = "Type I error (7 clean items)",
        col.main = "#1a1a19", xlab = "share of replications flagged")
abline(v = a, lty = 2, col = "#8a8a85")
legend("topright", c("50 + 50", "150 + 150 (glmer LRT not run)"), fill = cols, border = NA, bty = "n",
       cex = .8, title = "examinees per group", inset = c(0, .02))
dev.off()
cat("Wrote results/power/pilot_*.csv and results/figures/fig3_power_pilot.png\n")
