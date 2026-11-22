## The report's numbers come from the result files.
rep_md <- readLines(file.path("report", "report.md"), encoding = "UTF-8")
tab <- readLines(file.path("results", "report_tables.md"), encoding = "UTF-8")
rows <- tab[grepl("^\\| ", tab)]
check(sprintf("all %d generated table rows (Tables 3, 4, 8, 9, 12) appear verbatim in report.md", length(rows)),
      all(rows %in% rep_md))
r <- read.csv(file.path("results", "power_study", "rates.csv"))
m <- read.csv(file.path("results", "power_study", "rates_mcse.csv"))
check("rates_mcse.csv reproduces rates.csv (power and Type I to within .003)",
      max(abs(m$power - m$power_rep), na.rm = TRUE) < .003 && max(abs(m$type1 - m$type1_rep), na.rm = TRUE) < .003)
check("all 40 conditions, 300 replications each", length(unique(r$cell)) == 40 && all(r$reps == 300))
s1 <- read.csv(file.path("results", "sim_summary_rep1000_MIMIC-LRT-TREE.csv"))
check("Simulation 1: 1,000 replications per method, no failed fits", all(s1$n_rep == 1000) && all(s1$fits_na == 0))
fp <- read.csv(file.path("results", "power_study", "fit_problems.csv"))
check("Simulation 2: no failed fit in any method", all(fp$value[grepl("NA|fallback", fp$what)] == 0))
check("Figures 1-5 exist", all(file.exists(file.path("report", "figures", sprintf("figure%d.png", 1:5)))))
figs <- regmatches(rep_md, regexpr("figures/[A-Za-z0-9_]+\\.png", rep_md))
check("every figure the report cites exists", length(figs) == 5 && all(file.exists(file.path("report", figs))))
