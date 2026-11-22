## =============================================================================
## report_figures.R
## Figures 1-5 of the report, redrawn in one style from the result files in results/.
## No simulation is run. Base R graphics only.
##   Rscript R/report_figures.R
## Reads                                                     Writes
##   results/sim_summary_rep1000_MIMIC-LRT-TREE.csv          report/figures/figure1.png
##   results/sim_item_rates_rep1000_MIMIC-LRT-TREE.csv       report/figures/figure2.png
##   results/power/pilot_rates.csv, exact_items.csv,
##     exact_contamination_by_n.csv                          report/figures/figure3.png
##   results/power_study/rates.csv                           report/figures/figure4.png (Type I error)
##                                                           report/figures/figure5.png (power)
## Every method has one colour and one marker in every figure (validated categorical palette,
## markers as the second cue); benchmarks are grey. Titles are in the report's captions.
## =============================================================================

## The labels use UTF-8 characters (alpha, theta, mu, >=, en dash); switch the character type to
## UTF-8 if the session runs in a non-UTF-8 locale (e.g. "C" on a Linux server).
if (!isTRUE(l10n_info()[["UTF-8"]])) invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))
out <- file.path("report", "figures")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

## ---- Style ---------------------------------------------------------------------------------
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#8a8a85"; GRID <- "#e4e3df"; BAND <- "#f1f0ec"
COL <- c(mimic2011 = "#2a78d6", difLRT = "#eb6834", tree_gated = "#1baf7a", mh_purified = "#eda100",
         cml_purified = "#e87ba4", anchortest = "#008300", mimic_msp = "#6250d6", cml_robust5 = "#e34948")
PCH <- c(mimic2011 = 16, difLRT = 17, tree_gated = 15, mh_purified = 18,
         cml_purified = 1, anchortest = 0, mimic_msp = 2, cml_robust5 = 5)
LAB <- c(mimic2011 = "MIMIC, baseline specification", difLRT = "IRT-LRT (difLRT, purified)",
         tree_gated = "Rasch tree + anchortest", mh_purified = "Mantel–Haenszel, purified",
         cml_purified = "CML LRT, purified", anchortest = "anchortest without the gate",
         mimic_msp = "MIMIC with scale purification", cml_robust5 = "CML LRT, robust anchor (k = 5)")
f3 <- function(x) sub("^0\\.", ".", sprintf("%.3f", x))
f2 <- function(x) sub("^0\\.", ".", sprintf("%.2f", x))

save_png <- function(name, width, height, draw) {
  png(file.path(out, name), width = width, height = height, units = "in", res = 300,
      pointsize = 9, bg = "white", type = "cairo")
  on.exit(dev.off())
  par(family = "sans", col.axis = INK2, col.lab = INK, fg = INK2, las = 1, tcl = -0.25,
      mgp = c(2.2, 0.45, 0), cex.axis = 0.9)
  draw()
  invisible(file.path(out, name))
}
panel_label <- function(txt) mtext(txt, side = 3, line = 0.4, adj = 0, font = 2, col = INK, cex = 0.95)
hgrid <- function(at) abline(h = at, col = GRID, lwd = 0.8)
vgrid <- function(at) abline(v = at, col = GRID, lwd = 0.8)

## ---- Figure 1: Simulation 1, power and Type I error by method --------------------------------
s1 <- read.csv(file.path("results", "sim_summary_rep1000_MIMIC-LRT-TREE.csv"))
s1$key <- c("mimic2011", "difLRT", "tree_gated")[match(substr(s1$method, 1, 4), c("MIMI", "IRT-", "Rasc"))]
stopifnot(!anyNA(s1$key))
s1_lab <- c(mimic2011 = "MIMIC, baseline specification\n(z ≥ 1.96, one-sided)",
            difLRT = "IRT-LRT\n(difLRT, purified)", tree_gated = "Rasch tree\n+ anchored Wald test")
save_png("figure1.png", 6.5, 2.6, function() {
  par(mfrow = c(1, 2), mar = c(3.4, 0.8, 1.8, 1), oma = c(0, 10.5, 0, 0))
  y <- rev(seq_len(nrow(s1)))
  dot_panel <- function(est, lo, hi, xmax, label, ref = NULL) {
    plot(NA, xlim = c(0, xmax), ylim = c(0.5, nrow(s1) + 0.5), axes = FALSE, xlab = "Proportion of tests flagged", ylab = "")
    at <- pretty(c(0, xmax)); vgrid(at)
    if (!is.null(ref)) {
      abline(v = ref, col = INK2, lty = 2, lwd = 1)
      text(ref, nrow(s1) + 0.5, "α = .05", pos = 4, cex = 0.8, col = INK2, xpd = NA)
    }
    segments(lo, y, hi, y, col = COL[s1$key], lwd = 2)
    points(est, y, pch = PCH[s1$key], col = COL[s1$key], cex = 1.4)
    text(hi, y, f3(est), pos = 4, cex = 0.85, col = INK, xpd = NA)
    axis(1, at = at, labels = f2(at), col = GRID, col.ticks = INK2)
    panel_label(label)
  }
  dot_panel(s1$power, s1$power_lo, s1$power_hi, 0.5, "A. Power (i1, i5, i10)")
  axis(2, at = y, labels = s1_lab[s1$key], tick = FALSE, outer = TRUE, line = -0.6, col.axis = INK, cex.axis = 0.9)
  dot_panel(s1$type1, s1$type1_lo, s1$type1_hi, 0.2, "B. Type I error (7 clean items)", ref = 0.05)
})

## ---- Figure 2: Simulation 1, flag rate of each item by method --------------------------------
it <- read.csv(file.path("results", "sim_item_rates_rep1000_MIMIC-LRT-TREE.csv"))
it$dif <- as.logical(it$dif)
fig2_keys <- c(MIMIC = "mimic2011", LRT = "difLRT", TREE = "tree_gated")
save_png("figure2.png", 6.5, 3.3, function() {
  par(mar = c(4.6, 3.6, 2.6, 0.6))
  k <- nrow(it); off <- c(MIMIC = -0.22, LRT = 0, TREE = 0.22)
  ymax <- ceiling(max(it[, c("MIMIC_hi", "LRT_hi", "TREE_hi")]) * 10) / 10
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(0, ymax), axes = FALSE, xlab = "", ylab = "Proportion of replications flagged")
  for (i in which(it$dif)) rect(i - 0.5, 0, i + 0.5, ymax, col = BAND, border = NA)
  at <- pretty(c(0, ymax)); hgrid(at)
  abline(h = 0.05, col = INK2, lty = 2, lwd = 1)
  for (m in names(fig2_keys)) {
    x <- seq_len(k) + off[[m]]; key <- fig2_keys[[m]]
    segments(x, it[[paste0(m, "_lo")]], x, it[[paste0(m, "_hi")]], col = COL[[key]], lwd = 1.6)
    points(x, it[[m]], pch = PCH[[key]], col = COL[[key]], cex = 1.1)
  }
  axis(2, at = at, labels = f2(at), col = GRID, col.ticks = INK2)
  axis(1, at = seq_len(k), labels = it$item, tick = FALSE, line = -0.6, col.axis = INK)
  b_txt <- ifelse(it$dif, sprintf("b %+.1f → %+.1f", it$b_ref, it$b_foc), sprintf("b %+.1f", it$b_ref))
  mtext(b_txt, side = 1, at = seq_len(k), line = 1.2, cex = 0.7, col = INK2)
  mtext("Shaded: DIF items (the rate is power). Unshaded: clean items (Type I error). Dashed: .05. Bars: 95% CI.",
        side = 1, line = 2.9, adj = 0, cex = 0.75, col = INK2)
  legend("top", inset = c(0, -0.11), horiz = TRUE, xpd = NA, bty = "n",
         legend = c("MIMIC, baseline specification", "IRT-LRT (difLRT, purified)", "Rasch tree + anchored Wald test"),
         col = COL[fig2_keys], pch = PCH[fig2_keys], lwd = 1.6, cex = 0.9, text.col = INK, seg.len = 1.2)
})

## ---- Figure 3: the pilot, simulated against predicted ----------------------------------------
pr <- read.csv(file.path("results", "power", "pilot_rates.csv"))
ex <- read.csv(file.path("results", "power", "exact_items.csv")); exd <- ex[ex$dif > 0, ]
cn <- read.csv(file.path("results", "power", "exact_contamination_by_n.csv"))
## Predictions. Clean and all-other anchors: exact values at 50 and 150 per group. Theta known and
## the group proportions: the exact z at 50 + 50, scaled by sqrt(n / 50) (Var(delta-hat) ~ 1/n).
pw_scaled <- function(z, n) mean(pnorm(z * sqrt(n / 50) - qnorm(.975)) + pnorm(-z * sqrt(n / 50) - qnorm(.975)))
pred <- function(method, n, what) {
  cl <- cn[cn$n_per_group == n, ]
  switch(method,
    oracle_theta = if (what == "power") pw_scaled(exd$z_known_theta, n) else .05,
    marginal = , mimic2011 = if (what == "power") pw_scaled(exd$z_marginal, n) else .05,
    mimic_oracle = , lrt_oracle = , cml_oracle = , mh_oracle_nocc =
      if (what == "power") cl$mean_power_clean_anchor else .05,
    mimic_anchored = , cml_none = if (what == "power") cl$mean_power_all_anchor else cl$mean_type1_all_anchor,
    NA_real_)
}
show <- c(oracle_theta = "Benchmark: θ known", marginal = "Benchmark: group proportions",
          mimic2011 = "MIMIC, baseline specification", mimic_anchored = "MIMIC, f ~ x, all other items as anchor",
          mimic_oracle = "MIMIC, f ~ x, true clean anchor", lrt_difR = "difLRT, purified (baseline)",
          lrt_oracle = "difLRT statistic, true clean anchor", cml_none = "CML LRT, all other items as anchor",
          cml_purified = "CML LRT, purified", cml_oracle = "CML LRT, true clean anchor",
          mh_purified_nocc = "Mantel–Haenszel, purified", mh_oracle_nocc = "Mantel–Haenszel, true clean anchor",
          tree_gated = "Rasch tree + anchortest", anchortest = "anchortest without the gate")
N_COL <- c("50" = "#86b6ef", "150" = "#184f95"); N_PCH <- c("50" = 16, "150" = 15)
get_rate <- function(design, m, what) { v <- pr[[what]][pr$design == design & pr$method == m]; if (length(v)) v else NA_real_ }
save_png("figure3.png", 6.5, 4.4, function() {
  layout(matrix(1:2, 1), widths = c(1.45, 1))
  y <- rev(seq_along(show))
  pilot_panel <- function(what, xmax, label, ref, ylab) {
    plot(NA, xlim = c(0, xmax), ylim = c(0.5, length(show) + 0.5), axes = FALSE, xlab = "Proportion of tests flagged", ylab = "")
    at <- pretty(c(0, xmax)); vgrid(at)
    abline(v = ref, col = INK2, lty = 2, lwd = 1)
    for (i in seq_along(show)) {
      m <- names(show)[i]
      for (n in c(50, 150)) {
        p <- pred(m, n, what)
        if (!is.na(p)) segments(p, y[i] - 0.32, p, y[i] + 0.32, col = MUTED, lwd = 2)
      }
      sim <- c(get_rate("study", m, what), get_rate("study_n150", m, what))
      points(sim, rep(y[i], 2), pch = N_PCH, col = N_COL, cex = 1.15)
    }
    axis(1, at = at, labels = f2(at), col = GRID, col.ticks = INK2)
    if (ylab) axis(2, at = y, labels = show, tick = FALSE, line = -0.6, col.axis = INK, cex.axis = 0.85)
    panel_label(label)
  }
  par(mar = c(3.4, 15, 1.8, 0.8))
  pilot_panel("power", 1, "A. Power (i1, i5, i10)", 0.8, TRUE)
  par(mar = c(3.4, 0.8, 1.8, 0.8))
  pilot_panel("type1", 0.2, "B. Type I error (7 clean items)", 0.05, FALSE)
  legend("bottomright", c("50 + 50, simulated", "150 + 150, simulated", "predicted (exact)"),
         pch = c(N_PCH, NA), col = c(N_COL, MUTED), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
         bty = "o", box.col = GRID, bg = "white", cex = 0.85, text.col = INK, inset = 0.01)
})

## ---- Figures 4 and 5: Simulation 2, Type I error and power by condition ---------------------
rt <- read.csv(file.path("results", "power_study", "rates.csv"))
STUDY <- c("cml_purified", "mh_purified", "cml_robust5", "mimic2011", "mimic_msp", "tree_gated", "anchortest")   # LRT, MIMIC, tree families
NS <- c(50, 100, 150, 250)
study_panel <- function(z, what, ylim, title, ylab = "", band = NULL, ref = NULL, bench = FALSE) {
  plot(NA, xlim = c(35, 265), ylim = ylim, axes = FALSE, xlab = "", ylab = ylab)
  at <- pretty(ylim); hgrid(at)
  if (!is.null(band)) rect(35, band[1], 265, band[2], col = BAND, border = NA)
  if (!is.null(ref)) abline(h = ref, col = INK2, lty = 2, lwd = 1)
  if (bench) for (b in c("oracle_theta", "cml_oracle")) {
    s <- z[z$method == b, ]; s <- s[order(s$n), ]
    lines(s$n, s[[what]], col = MUTED, lwd = 1.8, lty = if (b == "oracle_theta") 3 else 4)
  }
  for (m in STUDY) {
    s <- z[z$method == m, ]; s <- s[order(s$n), ]
    v <- ifelse(s[[what]] > ylim[2], NA, s[[what]])           # off-scale values are not drawn ...
    lines(s$n, v, col = COL[[m]], lwd = 1.4)
    points(s$n, v, pch = PCH[[m]], col = COL[[m]], cex = if (PCH[[m]] == 18) 1.25 else 0.95)
    hi <- which(s[[what]] > ylim[2])                           # ... but marked at the top and printed
    if (length(hi)) {
      points(s$n[hi], rep(ylim[2], length(hi)), pch = 24, col = COL[[m]], bg = COL[[m]], cex = 0.8, xpd = NA)
      text(s$n[hi], ylim[2], f2(s[[what]][hi]), pos = 1, offset = 0.45, cex = 0.75, col = INK, xpd = NA)
    }
  }
  axis(1, at = NS, col = GRID, col.ticks = INK2, cex.axis = 0.85)
  axis(2, at = at, labels = f2(at), col = GRID, col.ticks = INK2, cex.axis = 0.85)
  mtext(title, side = 3, line = 0.3, adj = 0, font = 2, col = INK, cex = 0.8)
}
cond_title <- function(dr, im) sprintf("%s, %s", if (dr == "one-way") "One-way" else if (dr == "mixed") "Mixed" else "No DIF",
                                       if (im == 0) "μ = 0" else "μ = 0.5")
study_legend <- function(bench, ncol) {
  legend("top", legend = c(LAB[STUDY], if (bench) c("Benchmark: θ known", "Benchmark: CML LRT, true clean anchor")),
         col = c(COL[STUDY], if (bench) rep(MUTED, 2)), pch = c(PCH[STUDY], if (bench) c(NA, NA)),
         lty = c(rep(1, length(STUDY)), if (bench) c(3, 4)), lwd = c(rep(1.4, length(STUDY)), if (bench) rep(1.8, 2)),
         pt.cex = ifelse(c(PCH[STUDY], if (bench) c(1, 1)) == 18, 1.25, 1),
         bty = "n", cex = 0.9, text.col = INK, ncol = ncol, seg.len = 2.2, inset = 0.02)
}

## Figure 4: Type I error. Row 1: no-DIF conditions (and the legend); rows 2-3: DIF of 0.8 and 0.5.
t1_top <- max(rt$type1[rt$method %in% setdiff(STUDY, "mimic2011")], na.rm = TRUE)
t1_lim <- c(0, max(0.15, ceiling(t1_top * 20) / 20))
save_png("figure4.png", 6.5, 6.2, function() {
  layout(rbind(c(1, 2, 3, 3), 4:7, 8:11))
  par(mar = c(2.4, 2.6, 1.6, 0.4), oma = c(2, 2, 0, 0))
  for (im in c(0, 0.5)) study_panel(rt[rt$direction == "none" & rt$impact == im, ], "type1", t1_lim,
                                    cond_title("none", im), band = c(0.035, 0.065), ref = 0.05)
  par(mar = c(0, 0, 0, 0)); plot.new(); study_legend(FALSE, 1)
  par(mar = c(2.4, 2.6, 1.6, 0.4))
  for (dif in c(0.8, 0.5)) for (dr in c("one-way", "mixed")) for (im in c(0, 0.5)) {
    study_panel(rt[rt$direction == dr & rt$impact == im & rt$dif == dif, ], "type1", t1_lim,
                paste0(cond_title(dr, im), ", DIF ", dif), band = c(0.035, 0.065), ref = 0.05)
  }
  mtext("Examinees per group", side = 1, outer = TRUE, line = 0.6, col = INK, cex = 0.9)
  mtext("Type I error (two-sided α = .05)", side = 2, outer = TRUE, line = 0.6, col = INK, cex = 0.9, las = 0)
})

## Figure 5: power. Rows: DIF of 0.8 and 0.5 logits; columns: direction x impact; legend below.
save_png("figure5.png", 6.5, 4.9, function() {
  layout(rbind(1:4, 5:8, c(9, 9, 9, 9)), heights = c(1, 1, 0.42))
  par(mar = c(2.4, 2.6, 1.6, 0.4), oma = c(0, 2, 0, 0))
  for (dif in c(0.8, 0.5)) for (dr in c("one-way", "mixed")) for (im in c(0, 0.5)) {
    study_panel(rt[rt$direction == dr & rt$impact == im & rt$dif == dif, ], "power", c(0, 1),
                paste0(cond_title(dr, im), ", DIF ", dif), ref = 0.8, bench = TRUE)
  }
  par(mar = c(0, 0, 1.6, 0)); plot.new()
  mtext("Examinees per group", side = 3, line = 0.2, col = INK, cex = 0.9)
  study_legend(TRUE, 3)
  mtext("Power (two-sided α = .05)", side = 2, outer = TRUE, line = 0.6, col = INK, cex = 0.9, las = 0, adj = 0.72)
})

cat("Wrote", paste(file.path(out, sprintf("figure%d.png", 1:5)), collapse = ", "), "\n")
