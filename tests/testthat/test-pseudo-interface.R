pseudo_public_data <- function() {
  set.seed(1181)
  x <- sweep(matrix(rnorm(48 * 4), 48, 4), 2, c(1, 3, .7, 2), "*")
  a <- runif(48, 0, 4)
  r <- 3 * exp(-((a - 2 - .6 * tanh(x[, 1] + x[, 2])) / .7)^2)
  list(x = x, a = a, r = r)
}
pseudo_public_fit <- function(z, ...) {
  args <- list(x = z$x, a = z$a, r = z$r, method = "pseudo_direct", ndr = 2,
               B.initial = matrix(c(1, .3, 0, .2, 0, .2, 1, .4), 4, 2),
               maxitr = 0, keep.data = TRUE, lambda = c(.1, .3))
  updates <- list(...)
  args[names(updates)] <- updates
  do.call(orthoDr_pdose, args)
}
test_that("pseudo-direct metadata, prediction and retention agree", {
  z <- pseudo_public_data()
  for (intercept in c(FALSE, TRUE)) {
    f <- pseudo_public_fit(z, intercept = intercept)
    expect_equal(f$itr, 0)
    expect_identical(f$converge, f$gradient.norm <= 1e-6)
    expect_equal(f$fn, f$reward.scale * sqrt(f$objective.scaled), tolerance = 1e-14)
    expect_equal(f$dose.bandwidth.scale, sd(z$a) * f$bw / sqrt(2), tolerance = 1e-14)
    expect_equal(crossprod(f$B), diag(2), tolerance = 1e-14)
    expect_length(predict(f, z$x[FALSE, , drop = FALSE])$pred, 0)
    expect_equal(drop(predict(f, z$x)$pred), c(drop(predict(f, z$x[1:10, ])$pred),
      drop(predict(f, z$x[11:48, ])$pred)), tolerance = 2e-12)
    restored <- unserialize(serialize(f, NULL))
    expect_identical(predict(restored, z$x)$pred, predict(f, z$x)$pred)
    if (intercept) {
      restored$intercept <- NULL
      expect_error(predict(restored, z$x[1, , drop = FALSE]), "intercept")
    }
    small <- pseudo_public_fit(z, intercept = intercept, keep.data = FALSE)
    expect_false(any(c("x", "a", "r", "kernel", "grid.targets", "grid.scores") %in% names(small)))
    expect_true(is.finite(small$gradient.norm))
  }
})
test_that("pseudo-direct intercept follows dose origin and units", {
  z <- pseudo_public_data()
  f <- pseudo_public_fit(z, intercept = TRUE)
  moved <- z; moved$a <- z$a + 1000
  scaled <- z; scaled$a <- 7 * z$a
  fm <- pseudo_public_fit(moved, intercept = TRUE)
  fs <- pseudo_public_fit(scaled, intercept = TRUE)
  expect_equal(fm$lambda.init, f$lambda.init)
  expect_equal(fs$lambda.init, f$lambda.init)
  expect_equal(drop(predict(fm, z$x)$pred), 1000 + drop(predict(f, z$x)$pred), tolerance = 3e-11)
  expect_equal(drop(predict(fs, z$x)$pred), 7 * drop(predict(f, z$x)$pred), tolerance = 3e-11)
  z$a <- z$a - min(z$a) + .1
  constant <- pseudo_public_fit(z, intercept = TRUE, K = 1, lambda = c(.3, .1))
  expect_true(all(constant$W == 0))
  expect_identical(constant$intercept, .1)
  expect_identical(constant$lambda.init, .3)
})
test_that("pseudo controls and constant rewards fail clearly", {
  z <- pseudo_public_data()
  expect_error(pseudo_public_fit(z, lambda.reward = .1), "only used")
  expect_error(pseudo_public_fit(z, intercept.bounds = c(0, 4)), "only used")
  for (control in list(list(gamma = -1), list(gamma = 1.1), list(eta = 1), list(gtol = 0)))
    expect_error(pseudo_public_fit(z, control = control), "control")
  for (gamma in c(0, 1))
    expect_true(is.finite(pseudo_public_fit(z, control = list(gamma = gamma))$fn))
  z$r[] <- .1
  expect_error(pseudo_public_fit(z, B.initial = NULL), "Constant rewards")
})
test_that("genuine saved dose objects retain their predictions", {
  for (method in c("pseudo_direct", "direct", "direct_intercept")) {
    saved <- readRDS(test_path("fixtures", paste0("r-", method, ".rds")))
    expect_equal(drop(predict(saved$fit, saved$testx)$pred), saved$pred, tolerance = 2e-12)
  }
})
