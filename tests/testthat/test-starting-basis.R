test_that("nearly dependent supplied starts are refined before optimization", {
  set.seed(654)
  x <- matrix(rnorm(160), 40, 4)
  y <- 0.7 * x[, 1] - x[, 2]^2 + rnorm(40, sd = 0.25)
  v <- c(1, 2, 3, 4)
  start <- cbind(v, v + 2^-22 * c(0, 1, 0, 0))
  dependent <- cbind(v, v)

  for (method in c("sir", "save", "phd", "local", "seff")) {
    fit <- orthoDr_reg(x, y, method, ndr = 2, B.initial = start, maxitr = 5, ncore = 1)
    expect_gt(fit$itr, 0)
    expect_lt(norm(crossprod(fit$B) - diag(2), "F"), 1e-10)
  }
  expect_error(orthoDr_reg(x, y, ndr = 2, B.initial = dependent, ncore = 1), "full rank")

  time <- exp(y)
  event <- rep(c(TRUE, FALSE, TRUE), length.out = length(y))
  for (method in c("dn", "dm")) {
    fit <- orthoDr_surv(x, time, event, method, ndr = 2, B.initial = start, maxitr = 5, ncore = 1)
    expect_gt(fit$itr, 0)
    # The solver fits in standardized coordinates. Undo the returned unit conversion.
    working <- sweep(fit$B, 1L, apply(x, 2L, sd), "*")
    working <- sweep(working, 2L, sqrt(colSums(working^2)), "/")
    expect_lt(norm(crossprod(working) - diag(2), "F"), 1e-10)
  }
  expect_error(orthoDr_surv(x, time, event, B.initial = dependent, ncore = 1), "full rank")

  A <- diag(c(1, 4, 2, 3))
  fn <- function(B) -0.5 * sum(B * (A %*% B))
  grad <- function(B) -A %*% B
  expect_output(fit <- ortho_optim(start, fn, grad, maxitr = 5), "Gram-Schmidt")
  expect_gt(fit$itr, 0)
  expect_lt(norm(crossprod(fit$B) - diag(2), "F"), 1e-10)
  expect_error(capture.output(ortho_optim(dependent, fn, grad)), "full rank")
})
