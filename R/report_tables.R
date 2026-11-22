## =============================================================================
## report_tables.R
## Replication-level Monte Carlo SEs (MCSEs) for Simulation 2, and the report's results tables
## (Tables 3, 4, 8 and 9), built only from files already in results/. No simulation is run.
##   Rscript R/report_tables.R
## Writes results/power_study/rates_mcse.csv, results/report_tables.md (the tables as inserted
## into report/report.md) and results/report_tables_log.txt (check against rates.csv,
## sessionInfo()). Base R only.
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
meth <- c("oracle_theta", "cml_oracle", "mimic2011", "mimic_msp", "mh_purified", "cml_purified",
          "cml_robust5", "tree_gated", "anchortest")
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
lab <- c(oracle_theta = "Benchmark: θ known", cml_oracle = "Benchmark: CML LRT, true clean anchor",
         mimic2011 = "MIMIC, baseline specification", mimic_msp = "MIMIC with scale purification (M-SP)",
         mh_purified = "Mantel–Haenszel, purified", cml_purified = "CML LRT, purified",
         cml_robust5 = "CML LRT, robust anchor (k = 5)", tree_gated = "Rasch tree + anchortest",
         anchortest = "anchortest without the tree's gate")

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

cellv <- function(dir, dif, imp, v, meth) { d <- m[m$direction == dir & m$dif == dif & m$impact == imp & m$method == meth, ]
  d <- d[order(d$n), ]; paste(f3(d[[v]]), collapse = " / ") }
panels_t1 <- list(c("none", 0, 0, "No DIF, μ = 0"), c("none", 0, 0.5, "No DIF, μ = 0.5"),
                  c("one-way", 0.8, 0, "One-way 0.8, μ = 0"), c("one-way", 0.8, 0.5, "One-way 0.8, μ = 0.5"),
                  c("mixed", 0.8, 0.5, "Mixed 0.8, μ = 0.5"))
panels_pw <- list(c("one-way", 0.8, 0, "One-way 0.8, μ = 0"), c("one-way", 0.8, 0.5, "One-way 0.8, μ = 0.5"),
                  c("mixed", 0.8, 0.5, "Mixed 0.8, μ = 0.5"), c("one-way", 0.5, 0, "One-way 0.5, μ = 0"),
                  c("mixed", 0.5, 0.5, "Mixed 0.5, μ = 0.5"))
tab <- function(key, panels, v) {
  add("", key, paste0("| Method | ", paste(sapply(panels, `[`, 4), collapse = " | "), " |"),
      paste0("|---|", paste(rep("---", length(panels)), collapse = "|"), "|"))
  for (k in meth) add(paste0("| ", lab[k], " | ", paste(sapply(panels, function(p)
    cellv(p[1], as.numeric(p[2]), as.numeric(p[3]), v, k)), collapse = " | "), " |"))
}
tab("TABLE_TYPE1", panels_t1, "type1")
tab("TABLE_POWER", panels_pw, "power")
mx <- function(v, sub) round(max(m[[v]][sub], na.rm = TRUE), 3)
add("", "MCSE_NOTE", sprintf("largest MCSE: Type I %s (no-DIF cells %s); power %s",
    mx("type1_mcse", m$method %in% meth), mx("type1_mcse", m$method %in% meth & m$dif == 0),
    mx("power_mcse", m$method %in% meth)))
writeLines(md, file.path(out, "report_tables.md"))
