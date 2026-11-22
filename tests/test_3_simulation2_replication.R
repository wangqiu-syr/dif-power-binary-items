## Simulation 2: regenerate replication 1 of condition 1 (50 per group, one-way DIF of 0.5, no
## impact) with every method and compare with the stored statistics (rounded to 4 decimals).
if (!skip_if_missing(c("lavaan", "psychotree", "psychotools", "partykit"))) {
  source(file.path("R", "08_power_study.R"))
  tmp <- file.path(tempdir(), "study_test"); dir.create(tmp, showWarnings = FALSE)
  unlink(file.path(tmp, "*"))
  cell <- study_cells()[1, ]
  one_rep(cell, 1L, tmp)
  new <- readRDS(file.path(tmp, "cell01_0001.rds"))
  old <- readRDS(file.path("results", "power_study", "replications.rds"))
  old <- old[old$cell == 1 & old$rep == 1, ]
  m <- merge(new[, c("method", "item", "stat", "p")], data.frame(method = as.character(old$method), item = as.character(old$item),
             stat_old = old$stat, p_old = old$p), by = c("method", "item"))
  check("every method and item of the stored replication is regenerated", nrow(m) == nrow(old) && nrow(m) == nrow(new))
  d <- abs(m$stat - m$stat_old)
  check(sprintf("statistics within 0.001 of the stored ones (largest difference %.5f)", max(d, na.rm = TRUE)),
        max(d, na.rm = TRUE) < 0.001 && identical(is.na(m$stat), is.na(m$stat_old)))
  check("decisions at .05 identical", identical(m$p < .05, m$p_old < .05))
  info <- read.csv(file.path("results", "power_study", "replications_info.csv"))
  check("seed as documented: 20261119 + 1000 x cell + replication", new$seed[1] == 20261119 + 1000 * 1 + 1 &&
        info$seed[info$cell == 1 & info$rep == 1] == new$seed[1])
}
