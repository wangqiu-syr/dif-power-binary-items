## =============================================================================
## 09_power_study_summary.R
## Summaries of the full power study (R/08_power_study.R; docs/power_study.md).
##   Rscript R/09_power_study_summary.R [runs_dir]
## Writes results/power_study/: rates.csv (cell x method x flag rule), item_rates.csv,
## fit_problems.csv, replications.rds (every statistic), replications_info.csv; figures
## results/figures/fig4_study_power.png and
## fig5_study_type1.png.
## =============================================================================
source(file.path("R", "08_power_study.R"))          # study_cells(); the main block does not run when sourced
args <- commandArgs(trailingOnly = TRUE)
runs <- if (length(args) >= 1) args[1] else file.path("results", "power", "study_runs")
out <- file.path("results", "power_study"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
f <- list.files(runs, pattern = "^cell.*\\.rds$", full.names = TRUE)
x <- do.call(rbind, lapply(f, readRDS))
cells <- study_cells()
x <- merge(x, cells[, c("cell", "n", "direction", "impact", "dif")], by = "cell")
names(x)[names(x) == "dif"] <- "dif_size"
## Every statistic, compact (4.7 MB): one row per cell x replication x method x item; the design
## columns come from study_cells() by cell. Per-replication fields go to replications_info.csv.
comp <- x[, c("cell", "rep", "method", "item", "stat", "p")]
comp$stat <- round(comp$stat, 4); comp$p <- signif(comp$p, 4)
comp$method <- factor(comp$method); comp$item <- factor(comp$item)
saveRDS(comp, file.path(out, "replications.rds"), compress = "xz")
info <- x[!duplicated(paste(x$cell, x$rep)), c("cell", "rep", "seed", "tree_root_p", "msp_passes", "robust5_fallback", "secs")]
info$tree_root_p <- signif(info$tree_root_p, 4); info$secs <- round(info$secs, 2)
write.csv(info[order(info$cell, info$rep), ], file.path(out, "replications_info.csv"), row.names = FALSE)
a <- .05

## The 2011 tree: flagged only if the global test splits
tg <- x[x$method == "anchortest", ]; tg$method <- "tree_gated"
tg$p <- ifelse(tg$tree_root_p < a, tg$p, 1); tg$stat <- ifelse(tg$tree_root_p < a, tg$stat, 0)
x <- rbind(x, tg)
x$is_dif <- x$true_shift != 0
x$dir <- ifelse(x$item == "bundle", 1, sign(x$true_shift))

LABELS <- c(oracle_theta = "Benchmark: theta known", cml_oracle = "Benchmark: CML LRT, true clean anchor",
  mimic2011 = "MIMIC as in 2011 (no f ~ x)", mimic_msp = "MIMIC, f ~ x, scale purification",
  mh_purified = "Mantel-Haenszel, purified (no continuity correction)", cml_purified = "CML LRT, standard purification",
  cml_robust5 = "CML LRT, robust anchor (5 items)", cml_robust3 = "CML LRT, robust anchor (3 items)",
  tree_gated = "Rasch tree (gate) + anchortest, as in the rebuild", anchortest = "anchortest, no gate",
  difLRT = "difLRT as in 2011 (n = 50 only)", bundle = "Bundle test of i1 + i5 + i10")

## Flag rules. two-sided: p < .05. one-sided: z >= 1.645 in the hypothesised direction (one-way and
## no-DIF cells; the hypothesis is "easier for the focal group"). BH: Benjamini-Hochberg at .05
## across the 10 items of one method in one replication.
x$flag_two <- x$p < a
x$flag_one <- ifelse(x$direction == "mixed" | x$method == "difLRT", NA, x$stat >= qnorm(1 - a))
x$p_bh <- NA_real_
items_only <- x$item != "bundle"
key <- paste(x$cell, x$rep, x$method)[items_only]
x$p_bh[items_only] <- ave(x$p[items_only], key, FUN = function(p) p.adjust(p, "BH"))
x$flag_bh <- x$p_bh < a

rate <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
grp <- split(x, list(x$cell, x$method), drop = TRUE)
rates <- do.call(rbind, lapply(grp, function(z) {
  d <- z[z$is_dif, ]; c0 <- z[!z$is_dif, ]
  fdr <- if (z$item[1] == "bundle") NA_real_ else {
    per <- tapply(seq_len(nrow(z)), z$rep, function(i) { fl <- z$flag_bh[i]; if (!any(fl, na.rm = TRUE)) 0 else
      sum(fl & !z$is_dif[i], na.rm = TRUE) / sum(fl, na.rm = TRUE) }); mean(per) }
  data.frame(cell = z$cell[1], n = z$n[1], direction = z$direction[1], impact = z$impact[1], dif = z$dif_size[1],
             method = z$method[1], reps = length(unique(z$rep)),
             power = if (nrow(d)) rate(d$flag_two) else NA_real_, type1 = if (nrow(c0)) rate(c0$flag_two) else NA_real_,
             power_one_sided = if (nrow(d)) rate(d$flag_one) else NA_real_,
             type1_one_sided = if (nrow(c0)) rate(c0$flag_one) else NA_real_,
             power_bh = if (nrow(d)) rate(d$flag_bh) else NA_real_, type1_bh = if (nrow(c0)) rate(c0$flag_bh) else NA_real_,
             fdr_bh = fdr, na = mean(is.na(z$p)))
}))
## size-adjusted power: the critical p at which the method flags 5% of items in the no-DIF cell with
## the same n and impact
nul <- x[x$direction == "none" & x$item != "bundle", ]
qtab <- aggregate(p ~ n + impact + method, nul, function(p) quantile(p, a, na.rm = TRUE))
names(qtab)[4] <- "p_crit"
rates <- merge(rates, qtab, by = c("n", "impact", "method"), all.x = TRUE)
rates$power_size_adjusted <- NA_real_
for (i in which(!is.na(rates$p_crit) & rates$direction != "none" & rates$method != "tree_gated")) {
  z <- x[x$cell == rates$cell[i] & x$method == rates$method[i] & x$is_dif & x$item != "bundle", ]
  rates$power_size_adjusted[i] <- rate(z$p < rates$p_crit[i])
}
rates <- rates[order(rates$cell, match(rates$method, names(LABELS))), ]
rates$label <- unname(LABELS[rates$method])

item_rates <- aggregate(cbind(flag_two, flag_bh) ~ cell + method + item, x, mean, na.action = na.pass)
problems <- data.frame(
  what = c("MIMIC 2011 NA", "M-SP NA", "CML purified NA", "robust anchor fallback (k = 5)", "anchortest NA (excluding the first anchor item)",
           "M-SP passes: mean", "seconds per replication: mean"),
  value = c(rate(is.na(x$p[x$method == "mimic2011"])), rate(is.na(x$p[x$method == "mimic_msp"])),
            rate(is.na(x$p[x$method == "cml_purified"])), rate(tapply(x$robust5_fallback, paste(x$cell, x$rep), `[`, 1)),
            rate(tapply(is.na(x$p[x$method == "anchortest"]), paste(x$cell, x$rep)[x$method == "anchortest"], sum) > 1),
            rate(tapply(x$msp_passes, paste(x$cell, x$rep), `[`, 1)), rate(tapply(x$secs, paste(x$cell, x$rep), `[`, 1))))
rnd <- function(d) { i <- vapply(d, is.numeric, TRUE) & !names(d) %in% c("cell", "n", "reps"); d[i] <- lapply(d[i], round, 3); d }
write.csv(rnd(rates), file.path(out, "rates.csv"), row.names = FALSE)
write.csv(rnd(item_rates), file.path(out, "item_rates.csv"), row.names = FALSE)
write.csv(rnd(problems), file.path(out, "fit_problems.csv"), row.names = FALSE)
cat(sprintf("%d replications in %d cells\n", length(unique(paste(x$cell, x$rep))), length(unique(x$cell))))
print(rnd(problems), row.names = FALSE)

## ---- Figures ----------------------------------------------------------------------------------
## Seven methods in the reference categorical order (dataviz palette, validated), each with its own
## marker as secondary encoding; the theta-known benchmark in grey.
SHOW <- c(mimic2011 = "#2a78d6", mimic_msp = "#eb6834", mh_purified = "#1baf7a", cml_purified = "#eda100",
          cml_robust5 = "#e87ba4", tree_gated = "#008300", anchortest = "#6250d6")
PCH <- c(mimic2011 = 16, mimic_msp = 17, mh_purified = 15, cml_purified = 18, cml_robust5 = 1, tree_gated = 2, anchortest = 0)
SHORT <- c(mimic2011 = "MIMIC 2011", mimic_msp = "MIMIC, scale purification", mh_purified = "MH, purified",
           cml_purified = "CML LRT, purified", cml_robust5 = "CML LRT, robust anchor", tree_gated = "Rasch tree (gate)",
           anchortest = "anchortest, no gate")
## Rows: DIF of 0.8 and 0.5 logits; columns: direction x impact. Values above the y range (the 2011
## MIMIC's Type I error with impact) are printed at the top of the panel in the method's colour.
## Power panels add two grey references: theta known (dotted) and the CML LRT with the true clean
## anchor (dot-dash). Type I panels shade .035-.065 (the H1 criterion, .05 +/- .015).
panel_plot <- function(what, ylab, ylim, ref = NULL, band = NULL, file) {
  png(file, width = 2400, height = 1400, res = 170, type = "cairo")
  layout(matrix(c(1:4, 9, 5:8, 9), 2, byrow = TRUE), widths = c(1, 1, 1, 1, .95))
  par(mar = c(4.2, 4.3, 3, .6), las = 1, col.axis = "#3d3d3a", col.lab = "#3d3d3a", fg = "#8a8a85")
  for (dif in c(.8, .5)) for (dr in c("one-way", "mixed")) for (im in c(0, .5)) {
    z <- rates[rates$direction == dr & rates$impact == im & rates$dif == dif, ]
    plot(NA, xlim = c(40, 260), ylim = ylim, xaxt = "n", xlab = "examinees per group",
         ylab = if (dr == "one-way" && im == 0) sprintf("%s, DIF %.1f logits", ylab, dif) else "",
         main = sprintf("%s DIF of %.1f, %s", dr, dif, if (im == 0) "no impact" else "impact 0.5 SD"),
         col.main = "#1a1a19", cex.main = 1)
    axis(1, at = c(50, 100, 150, 250))
    if (!is.null(band)) rect(par("usr")[1], band[1], par("usr")[2], band[2], col = "#ececea", border = NA)
    if (!is.null(ref)) abline(h = ref, lty = 2, col = "#8a8a85")
    if (what == "power") for (bm in c("oracle_theta", "cml_oracle")) {
      b <- z[z$method == bm, ]; b <- b[order(b$n), ]
      lines(b$n, b[[what]], col = "#8a8a85", lwd = 2, lty = if (bm == "oracle_theta") 3 else 4) }
    for (m in names(SHOW)) { s <- z[z$method == m, ]; s <- s[order(s$n), ]
      lines(s$n, s[[what]], col = SHOW[[m]], lwd = 2); points(s$n, s[[what]], col = SHOW[[m]], pch = PCH[[m]], cex = 1.2)
      hi <- which(s[[what]] > ylim[2])
      if (length(hi)) text(s$n[hi], ylim[2], sprintf("%.2f", s[[what]][hi]), col = SHOW[[m]], pos = 1, cex = .85, font = 2) }
  }
  par(mar = c(4.2, 0, 3, 0)); plot.new()
  pw <- what == "power"
  legend("left", c(SHORT, if (pw) c("theta known (ceiling)", "CML LRT, clean anchor known")),
         col = c(SHOW, if (pw) rep("#8a8a85", 2)), pch = c(PCH, if (pw) c(NA, NA)), lwd = 2,
         lty = c(rep(1, length(SHOW)), if (pw) c(3, 4)), bty = "n", cex = .95)
  dev.off()
}
panel_plot("power", "power (i1, i5, i10)", c(0, 1), ref = .8, file = file.path("results", "figures", "fig4_study_power.png"))
t1max <- max(rates$type1[rates$method %in% setdiff(names(SHOW), "mimic2011") & rates$direction != "none"], na.rm = TRUE)
panel_plot("type1", "Type I error (7 clean items)", c(0, max(.15, ceiling(t1max * 20) / 20 + .02)), ref = a,
           band = c(.035, .065), file = file.path("results", "figures", "fig5_study_type1.png"))
cat("Wrote results/power_study/*.csv and results/figures/fig4_study_power.png, fig5_study_type1.png\n")
