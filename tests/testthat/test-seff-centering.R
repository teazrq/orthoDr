seff_centering_data <- function() {
  set.seed(7411)
  # Binary fractions make centering unchanged by these exact column shifts.
  x <- matrix(round(rnorm(64 * 4) * 16) / 16, 64, 4)
  x <- sweep(x, 2, colMeans(x), "-")
  x <- sweep(x, 2, c(1, 4, .5, 2), "*")
  shift <- c(16, -8, 4, 32)
  y <- 2 + x[, 1] + .4 * x[, 2]^2 + rnorm(64, sd = .2)
  B <- cbind(c(1, .5, -.25, 0), c(0, -.25, .5, 1))
  list(x = x, shifted = sweep(x, 2, shift, "+"), y = y, B = B)
}

test_that("SEFF is invariant to column shifts with default and supplied starts", {
  d <- seff_centering_data()
  fields <- c("B", "fn", "itr", "converge")
  for (ndr in 1:2) {
    for (supplied in c(FALSE, TRUE)) {
      B <- if (supplied) d$B[, seq_len(ndr), drop = FALSE] else NULL
      run <- function(x) orthoDr_reg(
        x, d$y, method = "seff", ndr = ndr, B.initial = B,
        bw = .7, maxitr = 5, ncore = 1
      )
      centered <- run(d$x)
      shifted <- run(d$shifted)
      expect_true(all(is.finite(shifted$B)))
      expect_true(is.finite(shifted$fn))
      expect_equal(shifted[fields], centered[fields], tolerance = 0)
    }
  }
})

test_that("SEFF uses the pooled scale and retains raw training data for prediction", {
  d <- seff_centering_data()
  bw <- .7
  B <- pracma::gramSchmidt(d$B)$Q
  fit <- orthoDr_reg(
    d$shifted, d$y, method = "seff", ndr = 2, B.initial = B,
    bw = bw, keep.data = TRUE, maxitr = 5, ncore = 1
  )
  Y <- as.matrix(scale(d$y) / nrow(d$x)^(-1 / 7) / sqrt(2))
  pooled <- sqrt(mean(apply(d$x, 2L, sd)^2))
  reference <- orthoDr:::seff_solver(
    pracma::gramSchmidt(B)$Q, d$x / pooled, Y, bw,
    1e-4, .2, .85, 1e-3, 1e-6, 1e-6, 1e-6, 1e-6,
    5L, 0L, 1L
  )
  reference$fn <- reference$fn * pooled^2
  fields <- c("B", "fn", "itr", "converge")
  expect_equal(fit[fields], reference[fields], tolerance = 2e-7)
  expect_identical(fit$x, d$shifted)
  expect_identical(fit$y, d$y)
  expect_identical(fit$bw, bw)

  testx <- d$shifted[c(2, 17, 41), , drop = FALSE] + .125
  denominator <- apply(d$shifted %*% fit$B, 2, sd) * bw * sqrt(2)
  expected <- vapply(seq_len(nrow(testx)), function(i) {
    delta <- sweep(d$shifted, 2, testx[i, ], "-") %*% fit$B
    weights <- exp(-rowSums(sweep(delta, 2, denominator, "/")^2))
    sum(weights * d$y) / sum(weights)
  }, numeric(1))
  expect_equal(predict(fit, testx)$pred, expected, tolerance = 1e-12)
  expect_identical(predict(unserialize(serialize(fit, NULL)), testx),
                   predict(fit, testx))
})
