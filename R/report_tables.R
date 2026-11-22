## =============================================================================
## report_tables.R
## Replication-level Monte Carlo SEs (MCSEs) for Simulation 2, and the report's results tables
## (Tables 3, 4, 8 and 9) and the methods-by-factors matrix (Table 12), built only from files
## already in results/. No simulation is run.
##   Rscript R/report_tables.R
## Writes results/power_study/rates_mcse.csv, results/power_study/factor_matrix.csv (the metrics
## and ratings behind Table 12), results/report_tables.md (the tables as inserted into
## report/report.md) and results/report_tables_log.txt (check against rates.csv, sessionInfo()).
## Base R only.
## =============================================================================
if (!isTRUE(l10n_info()[["UTF-8"]])) invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", "C.UTF-8")))
out <- "results"

rates <- read.csv("results/power_study/rates.csv")
rep   <- readRDS("results/power_study/replications.rds")
info  <- read.csv("results/power_study/replications_info.csv")
cells <- unique(rates[, c("cell", "n", "impact", "direction", "dif")])
dif_items <- c("i1", "i5", "i10")

# Rasch tree (gate) + anchortest: anchortest's flags count only when the tree's root test is
# significant (as in the power study)
at <- rep[rep$method == "anchortest", ]
at <- merge(at, info[, c("cell", "rep", "tree_root_p")], by = c("cell", "rep"))
at$p <- ifelse(at$tree_root_p < .05, at$p, 1)
at$method <- "tree_gated"
rep <- rbind(rep, at[, names(rep)])
rep <- rep[rep$item != "bundle", ]
rep$item <- as.character(rep$item)
rep <- merge(rep, cells[, c("cell", "dif")], by = "cell")
rep$is_dif <- rep$dif > 0 & rep$item %in% dif_items
rep <- rep[!is.na(rep$p), ]                 # untested (anchor) items are not in the denominator, as in rates.csv
rep$flag <- rep$p < .05

# per replication: number flagged (f) and number tested (t). The rate is the pooled ratio
# sum(f)/sum(t), as in rates.csv; its MCSE treats replications as clusters (ratio estimator):
# MCSE = sqrt(var(f - R t) / reps) / mean(t).
rep$one <- 1
agg <- aggregate(cbind(f = flag, t = one) ~ cell + method + rep + is_dif, data = rep, FUN = sum)
summ <- do.call(rbind, lapply(split(agg, list(agg$cell, agg$method, agg$is_dif), drop = TRUE),
  function(d) { R <- sum(d$f) / sum(d$t)
    data.frame(cell = d$cell[1], method = as.character(d$method[1]), is_dif = d$is_dif[1],
               rate = R, mcse = sqrt(var(d$f - R * d$t) / nrow(d)) / mean(d$t), reps = nrow(d)) }))
pw <- subset(summ, is_dif);  names(pw)[4:5] <- c("power_rep", "power_mcse")
t1 <- subset(summ, !is_dif); names(t1)[4:5] <- c("type1_rep", "type1_mcse")
m <- merge(rates, t1[, c("cell", "method", "type1_rep", "type1_mcse")], by = c("cell", "method"), all.x = TRUE)
m <- merge(m, pw[, c("cell", "method", "power_rep", "power_mcse")], by = c("cell", "method"), all.x = TRUE)
# bundle: one test per replication, binomial MCSE
b <- m$method == "bundle"
m$power_mcse[b] <- sqrt(m$power[b] * (1 - m$power[b]) / m$reps[b])
m$type1_mcse[b & m$dif == 0] <- sqrt(m$type1[b & m$dif == 0] * (1 - m$type1[b & m$dif == 0]) / m$reps[b & m$dif == 0])
m <- m[order(m$cell, m$method), ]
write.csv(m, file.path(out, "power_study", "rates_mcse.csv"), row.names = FALSE)

chk <- with(m[!b, ], c(max_abs_diff_type1 = max(abs(type1 - type1_rep), na.rm = TRUE),
                        max_abs_diff_power = max(abs(power - power_rep), na.rm = TRUE)))
cat("Check against rates.csv (rounded p in the .rds):\n"); print(round(chk, 4))
cat("Largest MCSE: Type I", round(max(m$type1_mcse, na.rm = TRUE), 4),
    "; power", round(max(m$power_mcse, na.rm = TRUE), 4), "\n")
cat("Largest MCSE, Type I in no-DIF cells:", round(max(m$type1_mcse[m$dif == 0 & !b], na.rm = TRUE), 4), "\n")
writeLines(capture.output(print(round(chk, 4)), sessionInfo()), file.path(out, "report_tables_log.txt"))

# compact tables for the report
## Methods in three families (LRT, MIMIC, Rasch tree), each starting with its baseline method; the
## theta-known benchmark first. difLRT ran only in the two n = 50, no-impact conditions.
blocks <- list(c(oracle_theta = ""),
               c(difLRT = "", cml_purified = "", mh_purified = "", cml_robust5 = "", cml_oracle = ""),
               c(mimic2011 = "", mimic_msp = ""),
               c(tree_gated = "", anchortest = ""))
names(blocks) <- c("", "LRT family", "MIMIC family", "Rasch tree family")
meth <- unlist(lapply(blocks, names), use.names = FALSE)
wide <- function(d, v) {
  d <- d[d$method %in% meth, ]
  w <- reshape(d[, c("method", "n", v)], idvar = "method", timevar = "n", direction = "wide")
  w[match(meth, w$method), ]
}
sel <- function(dir, dif, imp) m[m$direction == dir & m$dif == dif & m$impact == imp, ]
for (s in list(c("none", 0, 0), c("none", 0, 0.5), c("one-way", 0.8, 0), c("one-way", 0.8, 0.5),
               c("mixed", 0.8, 0), c("mixed", 0.8, 0.5), c("one-way", 0.5, 0), c("mixed", 0.5, 0))) {
  d <- sel(s[1], as.numeric(s[2]), as.numeric(s[3]))
  cat("\n==", s[1], "DIF", s[2], "impact", s[3], "| Type I (n=50,100,150,250)\n")
  print(wide(d, "type1"), row.names = FALSE)
  if (s[1] != "none") { cat("-- power\n"); print(wide(d, "power"), row.names = FALSE) }
}

# ---- Markdown tables for the report (values from the result files; written to report_tables.md)
f3 <- function(x) sub("^0\\.", ".", sprintf("%.3f", x))
f2 <- function(x) sub("^0\\.", ".", sprintf("%.2f", x))
md <- character(0)
add <- function(...) md <<- c(md, ...)
lab <- c(oracle_theta = "Benchmark: θ known", difLRT = "difLRT, purified (baseline; n = 50, μ = 0 only)",
         cml_purified = "CML LRT, purified", mh_purified = "Mantel–Haenszel, purified",
         cml_robust5 = "CML LRT, robust anchor (k = 5)", cml_oracle = "CML LRT, true clean anchor (benchmark)",
         mimic2011 = "MIMIC, baseline specification", mimic_msp = "MIMIC with scale purification (M-SP)",
         tree_gated = "Rasch tree + anchortest (baseline)", anchortest = "anchortest without the tree's gate")

s1 <- read.csv("results/sim_summary_rep1000_MIMIC-LRT-TREE.csv")
add("TABLE_SIM1", "| Method | Power (3 DIF items) | Type I error (7 clean items) |", "|---|---|---|")
s1_lab <- c("MIMIC, baseline specification (z \u2265 1.96, one-sided)", "IRT-LRT (difLRT, purified)", "Rasch tree + anchored Wald test")
s1_lab <- s1_lab[match(substr(s1$method, 1, 4), c("MIMI", "IRT-", "Rasc"))]
for (i in seq_len(nrow(s1))) add(sprintf("| %s | %s [%s, %s] | %s [%s, %s] |", s1_lab[i],
  f3(s1$power[i]), f3(s1$power_lo[i]), f3(s1$power_hi[i]), f3(s1$type1[i]), f3(s1$type1_lo[i]), f3(s1$type1_hi[i])))

ex <- read.csv("results/power/exact_items.csv"); d <- ex$dif > 0
ng <- read.csv("results/power/exact_n_for_power.csv")
add("", "TABLE_EXACT", "| | θ known | Item's own proportions (= baseline MIMIC) | Clean anchor | All other items as anchor |",
    "|---|---|---|---|---|",
    sprintf("| Power at 50 + 50, mean of i1, i5, i10 | %s | %s | %s | %s |", f2(mean(ex$power_known_theta[d])),
            f2(mean(ex$power_marginal[d])), f2(mean(ex$power_clean_anchor[d])), f2(mean(ex$power_all_anchor[d]))),
    sprintf("| Type I error, mean of the 7 clean items | .05 | .05 | .05 | %s |", f2(mean(ex$power_all_anchor[!d]))))
n80 <- ng[ng$table == "n for power .80" & ng$item %in% c("i1", "i5", "i10"), ]
rng <- function(v) paste0(min(v), "–", max(v))
add(sprintf("| n per group for .80 (range over i1, i5, i10) | %s | %s | %s | %s |", rng(n80$n_known_theta),
            rng(n80$n_marginal), rng(n80$n_clean_anchor), rng(n80$n_all_anchor)))

cellv <- function(dir, dif, imp, v, meth) {             # values at n = 50 / 100 / 150 / 250; "—" = not run
  d <- m[m$direction == dir & m$dif == dif & m$impact == imp & m$method == meth, ]
  x <- d[[v]][match(c(50, 100, 150, 250), d$n)]
  if (all(is.na(x))) return("—")
  paste(ifelse(is.na(x), "—", f3(x)), collapse = " / ") }
panels_t1 <- list(c("none", 0, 0, "No DIF, μ = 0"), c("none", 0, 0.5, "No DIF, μ = 0.5"),
                  c("one-way", 0.8, 0, "One-way 0.8, μ = 0"), c("one-way", 0.8, 0.5, "One-way 0.8, μ = 0.5"),
                  c("mixed", 0.8, 0.5, "Mixed 0.8, μ = 0.5"))
panels_pw <- list(c("one-way", 0.8, 0, "One-way 0.8, μ = 0"), c("one-way", 0.8, 0.5, "One-way 0.8, μ = 0.5"),
                  c("mixed", 0.8, 0.5, "Mixed 0.8, μ = 0.5"), c("one-way", 0.5, 0, "One-way 0.5, μ = 0"),
                  c("mixed", 0.5, 0.5, "Mixed 0.5, μ = 0.5"))
tab <- function(key, panels, v) {
  add("", key, paste0("| Method | ", paste(sapply(panels, `[`, 4), collapse = " | "), " |"),
      paste0("|---|", paste(rep("---", length(panels)), collapse = "|"), "|"))
  for (b in seq_along(blocks)) {
    if (nzchar(names(blocks)[b])) add(paste0("| **", names(blocks)[b], "** |", strrep(" |", length(panels))))
    for (k in names(blocks[[b]])) add(paste0("| ", lab[k], " | ", paste(sapply(panels, function(p)
      cellv(p[1], as.numeric(p[2]), as.numeric(p[3]), v, k)), collapse = " | "), " |"))
  }
}
tab("TABLE_TYPE1", panels_t1, "type1")
tab("TABLE_POWER", panels_pw, "power")
mx <- function(v, sub) round(max(m[[v]][sub], na.rm = TRUE), 3)
add("", "MCSE_NOTE", sprintf("largest MCSE: Type I %s (no-DIF cells %s); power %s",
    mx("type1_mcse", m$method %in% meth), mx("type1_mcse", m$method %in% meth & m$dif == 0),
    mx("power_mcse", m$method %in% meth)))
## ---- Table 12: nine methods x twelve factors (Yes / Semi / No) -----------------------------------
## Each factor of Table 5 has a rule (README of the report, Table 12 note). Statistic- and
## anchor-based factors (A2, A3, A4, B3, B4) are rated for the gated tree from its item tests without
## the gate (the anchortest rows), because the gate is a separate decision rule. difLRT is rated from
## the two n = 50, no-impact conditions in which it ran; "—" = not run.
fm <- setdiff(meth, "oracle_theta")
R <- rates
sa <- function(d) ifelse(is.na(d$power_size_adjusted), d$power, d$power_size_adjusted)
grade <- function(x, yes, semi) if (is.na(x)) "—" else if (x <= yes) "Yes" else if (x <= semi) "Semi" else "No"
bench <- R[R$method == "cml_oracle", ]
metric <- function(mm) {
  src <- function(f) if (mm == "tree_gated" && f %in% c("A2", "A3", "A4", "B3", "B4")) "anchortest" else mm
  get <- function(f) R[R$method == src(f), ]
  d <- get("A1"); a1 <- if (mm == "difLRT") NA else sum(d$dif == 0.8 & d$n == 250 & d$direction != "none" & d$power >= .8)
  d <- get("A2"); b <- bench[match(d$cell, bench$cell), ]; ow <- d$direction == "one-way"
  a2 <- mean((sa(b) - sa(d))[ow], na.rm = TRUE)
  d <- get("A3")
  a3 <- if (mm == "difLRT") { cp <- R[R$method == "cml_purified" & R$cell %in% d$cell & R$direction != "none", ]
          sa(cp) - sa(d[d$direction != "none", ]) } else {
          b <- bench[match(d$cell, bench$cell), ]; mean((sa(b) - sa(d))[d$direction == "mixed" & d$impact == 0], na.rm = TRUE) }
  d <- get("A5"); ow0 <- d$direction == "one-way" & d$impact == 0 & d$n <= 150
  a5g <- mean((d$power_one_sided - d$power)[ow0], na.rm = TRUE)
  a5t <- suppressWarnings(max(d$type1_one_sided[d$direction %in% c("none", "one-way")], na.rm = TRUE))
  if (!is.finite(a5t)) a5t <- NA                          # difLRT: unsigned statistic, no one-sided rule
  d <- get("B1"); b1 <- if (mm == "difLRT") NA else max(d$type1[d$n == 250], na.rm = TRUE)
  d <- get("B2"); b2 <- if (mm == "difLRT") NA else {
    n0 <- d[d$direction == "none" & d$impact == 0, ]; n5 <- d[d$direction == "none" & d$impact == 0.5, ]
    max(abs(n5$type1 - n0$type1[match(n5$n, n0$n)])) }
  d <- get("B3"); b3 <- max(d$type1[d$direction == "one-way" & d$impact == 0], na.rm = TRUE)
  d <- get("B4"); b4 <- max(d$type1[d$direction == "none" & d$impact == 0], na.rm = TRUE)
  d <- get("B6"); b6 <- max(d$fdr_bh[d$direction == "one-way"], na.rm = TRUE)
  data.frame(method = mm, A1_conditions_at_80 = a1, A2_shortfall_oneway = a2, A3_shortfall_mixed = a3,
             A5_gain_one_sided = a5g, A5_type1_one_sided_max = a5t, B1_type1_n250_max = b1,
             B2_impact_shift_max = b2, B3_type1_oneway_max = b3, B4_type1_nodif_max = b4, B6_fdr_oneway_max = b6)
}
fx <- do.call(rbind, lapply(fm, metric))
## Structural ratings: A4 (a joint test of the suspected items; shown with the CML statistic, the
## bundle test; possible but not run for the others, marked with a dagger), A6 and B5 (the gate).
a4 <- c(difLRT = "Semi†", cml_purified = "Yes", mh_purified = "Semi†", cml_robust5 = "Yes", cml_oracle = "Yes",
        mimic2011 = "Semi†", mimic_msp = "Semi†", tree_gated = "Semi†", anchortest = "Semi†")
rate <- data.frame(method = fx$method,
  A1 = sapply(fx$A1_conditions_at_80, function(k) if (is.na(k)) "—" else if (k == 4) "Yes" else if (k >= 1) "Semi" else "No"),
  A2 = sapply(fx$A2_shortfall_oneway, grade, .03, .12),
  A3 = sapply(fx$A3_shortfall_mixed, grade, .03, .08),
  A4 = unname(a4[fx$method]),
  A5 = ifelse(fx$method == "difLRT", "No", ifelse(fx$A5_gain_one_sided >= .05 & fx$A5_type1_one_sided_max <= .065, "Yes",
              ifelse(fx$A5_gain_one_sided >= .02, "Semi", "No"))),
  A6 = ifelse(fx$method == "tree_gated", "No", "Yes"),
  B1 = sapply(fx$B1_type1_n250_max, grade, .065, .10),
  B2 = sapply(fx$B2_impact_shift_max, grade, .02, .05),
  B3 = sapply(fx$B3_type1_oneway_max, grade, .065, .10),
  B4 = sapply(fx$B4_type1_nodif_max, grade, .065, .10),
  B5 = ifelse(fx$method == "tree_gated", "Yes", "No"),
  B6 = sapply(fx$B6_fdr_oneway_max, grade, .075, .15))
num <- sapply(fx, is.numeric); fx[num] <- lapply(fx[num], round, 3)
write.csv(merge(fx, rate, by = "method", sort = FALSE), file.path(out, "power_study", "factor_matrix.csv"), row.names = FALSE)
short <- c(difLRT = "difLRT, purified (baseline)ᵃ", cml_purified = "CML LRT, purified", mh_purified = "Mantel–Haenszel, purified",
           cml_robust5 = "CML LRT, robust anchor (k = 5)", cml_oracle = "CML LRT, true clean anchor (benchmark)",
           mimic2011 = "MIMIC, baseline specification", mimic_msp = "MIMIC with scale purification (M-SP)",
           tree_gated = "Rasch tree + anchortest (baseline)ᵇ", anchortest = "anchortest without the tree's gate")
fac <- c("A1", "A2", "A3", "A4", "A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6")
add("", "TABLE_MATRIX", paste0("| Method | ", paste(fac, collapse = " | "), " |"), paste0("|---|", strrep("---|", length(fac))))
for (b in 2:length(blocks)) {
  add(paste0("| **", names(blocks)[b], "** |", strrep(" |", length(fac))))
  for (k in names(blocks[[b]])) add(paste0("| ", short[k], " | ", paste(unlist(rate[rate$method == k, fac]), collapse = " | "), " |"))
}
writeLines(md, file.path(out, "report_tables.md"))
