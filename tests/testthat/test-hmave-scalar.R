# Independent reference using all-observation population variance for whitening
# and the sample variance of observed event times for the leave-one-out score.
hmave_scalar_cv <- function(x, y, event) {
  predictor <- as.numeric(x)
  event_x <- predictor[event]
  event_y <- y[event]
  h2 <- 2 * length(y)^(-1 / 2)
  predictor_variance <- mean((predictor - mean(predictor))^2)
  weights <- exp(-outer(event_x, event_x, "-")^2 / predictor_variance / h2)
  diag(weights) <- 0
  targets <- exp(-outer(event_y, event_y, "-")^2 / var(event_y) / h2)
  predicted <- sweep(weights %*% targets, 1L, rowSums(weights), "/")
  sum((targets - predicted)^2)
}

test_that("scalar hMave preserves whitening, event-only score and matrix shape", {
  set.seed(610217)
  n <- 36L
  x <- matrix(rnorm(n, sd = .35), ncol = 1L)
  y <- exp(.5 * x[, 1L] + rnorm(n, sd = .55))
  event <- rep(c(TRUE, TRUE, FALSE, TRUE), length.out = n)
  expected <- hmave_scalar_cv(x, y, event)

  expect_lt(mean((x - mean(x))^2), 1)
  large_variance <- mean((5 * x - mean(5 * x))^2)
  expect_gt(large_variance, 1)
  expect_gt(abs(large_variance - round(large_variance)), .01)

  for (multiplier in c(1, 2, 5, -5)) {
    transformed <- multiplier * x + 7
    reference <- hmave_scalar_cv(transformed, y, event)
    expect_equal(reference, expected, tolerance = 1e-11)
    for (start in list(NULL, matrix(-2, 1L, 1L), matrix(.6, 1L, 1L))) {
      fit <- hMave(transformed, y, event, 1L, B0 = start)
      expect_identical(dim(fit$B), c(1L, 1L))
      expect_equal(abs(unname(fit$B[1L, 1L])), 1, tolerance = 1e-12)
      expect_true(is.finite(fit$cv))
      expect_equal(fit$cv, reference, tolerance = 1e-11)
    }
  }

  expect_equal(hMave(x, as.matrix(y), as.matrix(as.numeric(event)), 1L),
               hMave(x, y, event, 1L), tolerance = 1e-12)
})

test_that("scalar hMave rejects undefined inputs before fitting", {
  x <- matrix(seq_len(12L), ncol = 1L)
  y <- exp(seq_len(12L) / 12)
  event <- rep(c(TRUE, FALSE), 6L)
  for (start in list(NULL, matrix(1, 1L, 1L))) {
    expect_error(hMave(matrix(4, 12L, 1L), y, event, 1L, start),
                 "positive finite variance")
    for (bad in c(NA_real_, Inf, NaN)) {
      expect_error(hMave(replace(x, 1L, bad), y, event, 1L, start),
                   "finite numeric values")
      expect_error(hMave(x, replace(y, 1L, bad), event, 1L, start),
                   "finite values")
    }
    expect_error(hMave(matrix(1, 1L, 1L), 1, TRUE, 1L, start),
                 "At least two observations")
  }
  expect_error(hMave(matrix(numeric(), 0L, 1L), numeric(), logical(), 1L),
               "At least two observations")
  expect_error(hMave(matrix("x", 12L, 1L), y, event, 1L), "finite numeric values")
  for (bad in list(NULL, 0, 2, .5, Inf, NA_real_, c(1, 1), "1"))
    expect_error(hMave(x, y, event, bad), "m0 must be 1")
  expect_error(hMave(x, y, event, 1L, B0 = 1), "B0 must be a matrix")
  expect_error(hMave(x, y, event, 1L, B0 = matrix(1, 1L, 2L)), "1 x 1 matrix")
  for (bad in c(0, Inf, NA_real_, NaN))
    expect_error(hMave(x, y, event, 1L, B0 = matrix(bad, 1L, 1L)),
                 "finite, nonzero")
  expect_error(hMave(x, y, event, 1L, B0 = matrix("1", 1L, 1L)),
               "finite, nonzero")
})

test_that("hMave requires at least two distinct observed event times at every p", {
  y <- seq_len(12L)
  for (p in c(1L, 2L)) {
    x <- cbind(y, sin(y))[, seq_len(p), drop = FALSE]
    for (event in list(rep(FALSE, 12L), seq_len(12L) == 1L,
                       seq_len(12L) <= 2L)) {
      times <- y
      if (sum(event) == 2L) times[event] <- 1
      expect_error(hMave(x, times, event, 1L),
                   "At least two distinct observed event times")
    }
  }
})
