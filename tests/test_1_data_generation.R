## Data generation: layout, determinism, and the original order of random draws.
source(file.path("R", "functions.R"))
b_ref <- c(-1, -1, -1, 0, 0, 0, 0, 1, 1, 1)
b_foc <- c(-1.8, -1, -1, 0, -0.8, 0, 0, 1, 1, 0.2)
set.seed(101); d1 <- simulate_dif_data(b_ref, b_foc)
set.seed(101); d2 <- simulate_dif_data(b_ref, b_foc)
check("50 + 50 rows, group column x then i1-i10", identical(dim(d1), c(100L, 11L)) && identical(names(d1), c("x", ITEMS)))
check("group coded 0 (reference, first 50) and 1 (focal)", identical(d1$x, rep(0:1, each = 50)))
check("responses are 0/1", all(unlist(d1[, ITEMS]) %in% 0:1))
check("the same seed gives the same data set", identical(d1, d2))
## the original order: reference abilities, reference responses (one runif vector, persons within
## items), then the focal group
set.seed(101)
th_r <- rnorm(50); u_r <- runif(500); th_f <- rnorm(50); u_f <- runif(500)
ref <- matrix(as.integer(u_r < plogis(outer(th_r, b_ref, "-"))), 50)
foc <- matrix(as.integer(u_f < plogis(outer(th_f, b_foc, "-"))), 50)
check("draw order matches the original respon.1pl() code", all(as.matrix(d1[, ITEMS]) == rbind(ref, foc)))
## Simulation 2's generator: no impact and no DIF gives the same draws as Simulation 1's
source(file.path("R", "08_power_study.R"))
cl <- study_cells(); cell <- cl[cl$direction == "none" & cl$impact == 0 & cl$n == 50, ]
set.seed(7); a <- simulate_cell(cell)$dat; set.seed(7); b <- simulate_dif_data(b_ref, b_ref)
check("Simulation 2's generator equals Simulation 1's when there is no DIF and no impact", identical(a, b))
check("40 conditions, 8 of them without DIF; difLRT in 2", nrow(cl) == 40 && sum(cl$direction == "none") == 8 && sum(cl$diflrt) == 2)
