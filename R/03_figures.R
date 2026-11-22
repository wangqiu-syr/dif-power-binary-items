## =============================================================================
## 03_figures.R
## Merges the parallel simulation runs and draws the figures: MIMIC vs. IRT-LRT
## vs. Rasch tree.
##
## Step 1, merge. If per-seed files results/sim_replications_rep*_seed*_MIMIC-LRT-TREE.csv
## exist (e.g. 8 runs of R/02_simulation.R 125 <seed>, seeds 101-108), they are
## combined into one set of files for all replications:
##   results/sim_replications_rep1000_MIMIC-LRT-TREE.csv  every statistic and flag, with a seed column
##   results/sim_item_rates_rep1000_MIMIC-LRT-TREE.csv    flag rate per item and method, 95% CI
##   results/sim_summary_rep1000_MIMIC-LRT-TREE.csv       power and Type I error, 95% CI
## (1000 = the total number of replications). The per-seed files can then be
## deleted; without them the script reads the merged file instead.
##
## Step 2, figures (base R graphics only; no extra packages):
##   results/figures/fig1_power_type1.{png,pdf}      power and Type I error by method
##   results/figures/fig2_item_flag_rates.{png,pdf}  flag rate per item by method
##
## Run from the repository root:
##   Rscript R/03_figures.R
## or from the R / RStudio console:  source("R/03_figures.R")
## =============================================================================

TAG <- "MIMIC-LRT-TREE"
seed_files <- list.files("results", pattern = sprintf("^sim_replications_rep[0-9]+_seed[0-9]+_%s\\.csv$", TAG),
                         full.names = TRUE)
if (length(seed_files)) {
  ## Read every seed; the seed column keeps replication ids unique.
  long <- do.call(rbind, lapply(seed_files, function(f) {
    d <- read.csv(f)
    cbind(seed = as.integer(sub(".*_seed([0-9]+)_.*", "\\1", f)), d)
  }))
  ## Item parameters (identical in every run) come from a per-seed item-rates file.
  items <- read.csv(sub("sim_replications_", "sim_item_rates_", seed_files[1]))[, c("item", "b_ref", "b_foc", "dif")]
  n_rep <- nrow(unique(long[, c("seed", "rep")]))
  write.csv(long, file.path("results", sprintf("sim_replications_rep%d_%s.csv", n_rep, TAG)), row.names = FALSE)
  cat(sprintf("Merged %d per-seed files into sim_replications_rep%d_%s.csv\n", length(seed_files), n_rep, TAG))
} else {
  merged <- list.files("results", pattern = sprintf("^sim_replications_rep[0-9]+_%s\\.csv$", TAG), full.names = TRUE)
  if (length(merged) != 1) stop("Expected per-seed files or one merged sim_replications_rep*_", TAG, ".csv in results/")
  long  <- read.csv(merged)
  items <- read.csv(sub("sim_replications_", "sim_item_rates_", merged))[, c("item", "b_ref", "b_foc", "dif")]
  n_rep <- nrow(unique(long[, c("seed", "rep")]))
}
long$rid <- paste(long$seed, long$rep)
dif_items <- items$item[items$dif]
long$dif <- long$item %in% dif_items

METHODS <- c("MIMIC", "LRT", "TREE")
LABELS  <- c(MIMIC = "MIMIC (z >= 1.96)", LRT = "IRT-LRT (difLRT, purified)",
             TREE = "Rasch tree + anchored Wald")
COLS    <- c(MIMIC = "#2a78d6", LRT = "#eb6834", TREE = "#1baf7a")   # validated categorical slots 1-3
PCH     <- c(MIMIC = 16, LRT = 17, TREE = 15)                        # shape as a second cue
INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e4e3df"; BAND <- "#f1f0ec"

cat(sprintf("%d replications (seeds %s)\n", n_rep, paste(sort(unique(long$seed)), collapse = ", ")))

## ---- Pooled power / Type I error ------------------------------------------------------
## The replication is the independent unit, so the CI uses the spread of the
## per-replication rates (failed fits, flag = NA, are left out as in 02).
rate_ci <- function(d) {
  per_rep <- tapply(d$flag, d$rid, mean, na.rm = TRUE)
  per_rep <- per_rep[is.finite(per_rep)]
  m <- mean(per_rep); se <- sd(per_rep) / sqrt(length(per_rep))
  c(est = m, lo = max(0, m - 1.96 * se), hi = min(1, m + 1.96 * se), n_rep = length(per_rep))
}
pooled <- do.call(rbind, lapply(METHODS, function(mt) {
  d <- long[long$method == mt, ]
  p <- rate_ci(d[d$dif, ]); t1 <- rate_ci(d[!d$dif, ])
  data.frame(method = LABELS[[mt]], power = p[["est"]], power_lo = p[["lo"]], power_hi = p[["hi"]],
             type1 = t1[["est"]], type1_lo = t1[["lo"]], type1_hi = t1[["hi"]],
             n_rep = p[["n_rep"]], fits_na = sum(is.na(d$flag)))
}))
num <- sapply(pooled, is.numeric) & !names(pooled) %in% c("n_rep", "fits_na")
pooled[num] <- round(pooled[num], 4)
print(pooled, row.names = FALSE)
write.csv(pooled, file.path("results", sprintf("sim_summary_rep%d_%s.csv", n_rep, TAG)), row.names = FALSE)

## ---- Per-item flag rates (Wilson 95% CI over replications) ---------------------------------
wilson <- function(x) {
  x <- x[!is.na(x)]; n <- length(x); p <- mean(x); z <- 1.96
  ctr <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(est = p, lo = ctr - hw, hi = ctr + hw)
}
item_tab <- do.call(rbind, lapply(METHODS, function(mt) do.call(rbind, lapply(items$item, function(it) {
  w <- wilson(long$flag[long$method == mt & long$item == it])
  data.frame(method = mt, item = it, est = w[["est"]], lo = w[["lo"]], hi = w[["hi"]])
}))))

## Merged item-rates file: item parameters, then rate / lo / hi per method.
item_rates <- items
for (mt in METHODS) {
  d <- item_tab[item_tab$method == mt, ]
  item_rates[[mt]] <- round(d$est, 4)
  item_rates[[paste0(mt, "_lo")]] <- round(d$lo, 4)
  item_rates[[paste0(mt, "_hi")]] <- round(d$hi, 4)
}
write.csv(item_rates, file.path("results", sprintf("sim_item_rates_rep%d_%s.csv", n_rep, TAG)), row.names = FALSE)

## ---- Drawing helpers ---------------------------------------------------------------------
dir.create(file.path("results", "figures"), showWarnings = FALSE, recursive = TRUE)
save_both <- function(name, width, height, draw) {
  png(file.path("results", "figures", paste0(name, ".png")), width = width, height = height,
      units = "in", res = 200, bg = "white")
  draw(); dev.off()
  pdf(file.path("results", "figures", paste0(name, ".pdf")), width = width, height = height,
      bg = "white")
  draw(); dev.off()
}
base_par <- function(...) par(family = "sans", col.axis = INK2, col.lab = INK, fg = INK2,
                              las = 1, tcl = -0.25, mgp = c(2.4, 0.5, 0), ...)

## ---- Figure 1: power and Type I error by method --------------------------------------------
fig1 <- function() {
  base_par(mfrow = c(1, 2), mar = c(4, 1, 3, 1.2), oma = c(0, 13, 2.2, 0))
  y <- rev(seq_along(METHODS))
  panel <- function(est, lo, hi, title, xmax, ref = NULL) {
    plot(NA, xlim = c(0, xmax), ylim = c(0.5, length(METHODS) + 0.5), axes = FALSE,
         xlab = "Proportion of items flagged", ylab = "")
    abline(v = pretty(c(0, xmax)), col = GRID, lwd = 1)
    if (!is.null(ref)) {
      abline(v = ref, col = INK2, lty = 2, lwd = 1.2)
      text(ref, length(METHODS) + 0.45, "nominal .05", pos = 4, cex = 0.75, col = INK2, xpd = NA)
    }
    segments(lo, y, hi, y, col = COLS[METHODS], lwd = 2)
    points(est, y, pch = PCH[METHODS], col = COLS[METHODS], cex = 1.5)
    text(est, y + 0.28, sprintf("%.3f", est), cex = 0.8, col = INK, xpd = NA)
    axis(1, at = pretty(c(0, xmax)), col = GRID, col.ticks = INK2)
    mtext(title, side = 3, line = 0.8, adj = 0, font = 2, col = INK, cex = 0.95)
  }
  panel(pooled$power, pooled$power_lo, pooled$power_hi, "Power (3 DIF items)",
        max(0.5, ceiling(max(pooled$power_hi) * 10) / 10))
  axis(2, at = y, labels = LABELS[METHODS], tick = FALSE, outer = TRUE, line = -0.5, col.axis = INK, cex.axis = 0.9)
  panel(pooled$type1, pooled$type1_lo, pooled$type1_hi, "Type I error (7 non-DIF items)",
        max(0.2, ceiling(max(pooled$type1_hi) * 20) / 20), ref = 0.05)
  mtext(sprintf("MIMIC vs. IRT-LRT vs. Rasch tree: %d replications, N = 50 + 50", n_rep),
        side = 3, outer = TRUE, line = 0.6, adj = 0, font = 2, col = INK, cex = 1.05)
}
save_both("fig1_power_type1", 9, 3.6, fig1)

## ---- Figure 2: flag rate per item by method ------------------------------------------------
fig2 <- function() {
  base_par(mar = c(5.2, 4, 3.6, 1))
  k <- nrow(items); off <- c(MIMIC = -0.22, LRT = 0, TREE = 0.22)
  ymax <- ceiling(max(item_tab$hi) * 10) / 10
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(0, ymax), axes = FALSE, xlab = "", ylab = "Flag rate")
  for (i in which(items$dif)) rect(i - 0.5, 0, i + 0.5, ymax, col = BAND, border = NA)
  abline(h = pretty(c(0, ymax)), col = GRID, lwd = 1)
  abline(h = 0.05, col = INK2, lty = 2, lwd = 1.2)
  for (mt in METHODS) {
    d <- item_tab[item_tab$method == mt, ]; x <- seq_len(k) + off[[mt]]
    segments(x, d$lo, x, d$hi, col = COLS[[mt]], lwd = 2)
    points(x, d$est, pch = PCH[[mt]], col = COLS[[mt]], cex = 1.3)
  }
  axis(2, at = pretty(c(0, ymax)), col = GRID, col.ticks = INK2)
  axis(1, at = seq_len(k), labels = items$item, tick = FALSE, line = -0.4, col.axis = INK)
  shift <- ifelse(items$dif, sprintf("b %+.1f -> %+.1f", items$b_ref, items$b_foc), sprintf("b %+.1f", items$b_ref))
  mtext(shift, side = 1, at = seq_len(k), line = 1.4, cex = 0.7, col = INK2)
  mtext("Shaded = DIF item (rate is power); unshaded = non-DIF (rate is Type I error); dashed = .05",
        side = 1, line = 3.2, adj = 0, cex = 0.75, col = INK2)
  legend("top", inset = c(0, -0.13), horiz = TRUE, xpd = NA, bty = "n", legend = LABELS[METHODS],
         col = COLS[METHODS], pch = PCH[METHODS], lwd = 2, cex = 0.85, text.col = INK, seg.len = 1.2)
  mtext(sprintf("Flag rate per item, %d replications (95%% CI)", n_rep), side = 3, line = 2.2,
        adj = 0, font = 2, col = INK, cex = 1.05)
}
save_both("fig2_item_flag_rates", 9, 4.6, fig2)

cat("Wrote results/figures/fig1_power_type1.{png,pdf}, results/figures/fig2_item_flag_rates.{png,pdf},\n",
    sprintf("      results/sim_{item_rates,summary}_rep%d_%s.csv\n", n_rep, TAG), sep = "")
