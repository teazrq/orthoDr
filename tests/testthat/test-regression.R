# Tests for orthoDr regression: orthoDr_reg, hMave, pSAVE, and predict
# Covers all 5 methods, parameters, edge cases, and prediction

# ── Helpers ───────────────────────────────────────────────────────────────────

make_regression_data <- function(n = 80, p = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  testx <- matrix(rnorm(20 * p), nrow = 20, ncol = p)
  list(x = x, testx = testx)
}

# ── 1. All 5 methods: basic fitting ──────────────────────────────────────────

test_that("orthoDr_reg fits all 5 methods and returns correct structure", {
  set.seed(123)
  d <- make_regression_data(n = 80, p = 4)

  y_sir   <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(80, sd = 0.2)
  y_save  <- d$x[, 1] - 0.5 * d$x[, 2] + rnorm(80, sd = 0.3)
  y_phd   <- d$x[, 1]^2 - 0.6 * d$x[, 2]^2 + rnorm(80, sd = 0.2)
  y_local <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(80, sd = 0.3)
  y_seff  <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(80, sd = 0.3)

  methods <- c("sir", "save", "phd", "local", "seff")
  responses <- list(y_sir, y_save, y_phd, y_local, y_seff)

  for (i in seq_along(methods)) {
    fit <- orthoDr_reg(
      d$x, responses[[i]], method = methods[[i]],
      ndr = 1, keep.data = TRUE, maxitr = 40, verbose = FALSE
    )

    # Class
    expect_equal(class(fit), c("orthoDr", "fit", "reg"))

    # Dimensions
    expect_equal(dim(fit$B), c(ncol(d$x), 1L))
    expect_true(all(is.finite(fit$B)))

    # Required components
    expect_true("fn"        %in% names(fit))
    expect_true("itr"       %in% names(fit))
    expect_true("converge"  %in% names(fit))
    expect_true("method"    %in% names(fit))
    expect_true("keep.data" %in% names(fit))
    expect_equal(fit$method,    methods[[i]])
    expect_equal(fit$keep.data, TRUE)

    # Prediction
    pred <- predict(fit, d$testx)
    expect_equal(class(pred), c("orthoDr", "predict", "reg"))
    expect_equal(length(pred$pred), nrow(d$testx))
    expect_true(all(is.finite(pred$pred)))
    expect_true(is.numeric(pred$pred))
  }
})

# ── 2. ndr > 1 ────────────────────────────────────────────────────────────────

test_that("orthoDr_reg works with ndr = 2", {
  set.seed(234)
  d <- make_regression_data(n = 100, p = 5)

  # A model driven by 2 directions
  y <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + 0.3 * d$x[, 3] + rnorm(100, sd = 0.3)

  for (m in c("sir", "phd", "save", "local", "seff")) {
    fit <- orthoDr_reg(
      d$x, y, method = m,
      ndr = 2, keep.data = TRUE, maxitr = 50, verbose = FALSE
    )
    expect_equal(dim(fit$B), c(ncol(d$x), 2L))
    expect_equal(class(fit), c("orthoDr", "fit", "reg"))
    expect_true(all(is.finite(fit$B)))

    # Prediction with ndr=2 should still work
    pred <- predict(fit, d$testx)
    expect_equal(length(pred$pred), nrow(d$testx))
    expect_true(all(is.finite(pred$pred)))
  }
})

test_that("orthoDr_reg works with ndr = 3", {
  set.seed(345)
  d <- make_regression_data(n = 120, p = 6)
  y <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(120, sd = 0.3)

  fit <- orthoDr_reg(
    d$x, y, method = "sir",
    ndr = 3, keep.data = TRUE, maxitr = 40, verbose = FALSE
  )
  expect_equal(dim(fit$B), c(ncol(d$x), 3L))
  pred <- predict(fit, d$testx)
  expect_equal(length(pred$pred), nrow(d$testx))
})

# ── 3. B.initial ──────────────────────────────────────────────────────────────

test_that("orthoDr_reg accepts custom B.initial", {
  set.seed(456)
  d <- make_regression_data(n = 80, p = 4)
  y <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(80, sd = 0.2)

  # Provide a valid B.initial (orthogonal)
  B0 <- pracma::gramSchmidt(matrix(rnorm(4 * 1), 4, 1))$Q
  fit <- orthoDr_reg(
    d$x, y, method = "sir",
    ndr = 1, B.initial = B0, maxitr = 40, verbose = FALSE
  )
  expect_equal(class(fit), c("orthoDr", "fit", "reg"))
  expect_equal(dim(fit$B), c(4L, 1L))
  expect_true(fit$converge)
})

test_that("orthoDr_reg rejects bad B.initial dimensions", {
  set.seed(567)
  d <- make_regression_data(n = 60, p = 4)
  y <- 1 + d$x[, 1] + rnorm(60)

  expect_error(
    orthoDr_reg(d$x, y, method = "sir", ndr = 2, B.initial = matrix(1:8, 4, 2)),
    NA  # should succeed
  )

  expect_error(
    orthoDr_reg(d$x, y, method = "sir", ndr = 1, B.initial = matrix(1:6, 3, 2)),
    "Dimension of B.initial is not correct"
  )
})

test_that("orthoDr_reg handles non-orthogonal B.initial gracefully", {
  set.seed(678)
  d <- make_regression_data(n = 80, p = 4)
  y <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(80, sd = 0.3)

  # Non-orthogonal B0 — orthoDr_reg calls gramSchmidt silently
  B0 <- matrix(rnorm(4 * 2), 4, 2)
  fit <- orthoDr_reg(
    d$x, y, method = "sir", ndr = 2,
    B.initial = B0, maxitr = 40, verbose = FALSE
  )
  expect_equal(class(fit), c("orthoDr", "fit", "reg"))
  expect_equal(dim(fit$B), c(4L, 2L))
  # Gram-Schmidt correction: B should be orthonormal despite non-orthogonal input
  BtB <- crossprod(fit$B)
  expect_equal(diag(as.matrix(BtB)), rep(1, 2), tolerance = 1e-5)
})

# ── 4. bw (bandwidth) ─────────────────────────────────────────────────────────

test_that("orthoDr_reg accepts custom bw", {
  set.seed(789)
  d <- make_regression_data(n = 80, p = 4)
  y <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(80, sd = 0.2)

  fit_default <- orthoDr_reg(
    d$x, y, method = "sir", ndr = 1, maxitr = 30, verbose = FALSE
  )
  fit_custom <- orthoDr_reg(
    d$x, y, method = "sir", ndr = 1, bw = 1.5, maxitr = 30, verbose = FALSE
  )

  expect_equal(class(fit_default), c("orthoDr", "fit", "reg"))
  expect_equal(class(fit_custom), c("orthoDr", "fit", "reg"))
  # Different bw may give different results — at minimum both converge
})

# ── 5. keep.data and prediction ──────────────────────────────────────────────

test_that("regression prediction requires keep.data = TRUE", {
  set.seed(890)
  d <- make_regression_data(n = 60, p = 3)
  y <- d$x[, 1]^2 + rnorm(60)

  fit <- orthoDr_reg(
    d$x, y, method = "phd",
    ndr = 1, keep.data = FALSE, maxitr = 30, verbose = FALSE
  )
  expect_error(
    predict(fit, d$testx),
    "Need the original data for prediction"
  )
})

test_that("predict validates testx input", {
  set.seed(901)
  d <- make_regression_data(n = 60, p = 3)
  y <- 1 + d$x[, 1] + rnorm(60)

  fit <- orthoDr_reg(
    d$x, y, method = "sir", ndr = 1, keep.data = TRUE,
    maxitr = 30, verbose = FALSE
  )

  expect_error(predict(fit), "testx is missing")
  expect_error(
    predict(fit, matrix(letters[1:9], 3, 3)),
    "testx must be a numerical matrix"
  )
  expect_error(
    predict(fit, as.data.frame(d$testx)),
    "testx must be a numerical matrix"
  )
})

test_that("prediction is consistent for same inputs", {
  set.seed(1012)
  d <- make_regression_data(n = 80, p = 4)
  y <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(80, sd = 0.2)

  fit <- orthoDr_reg(
    d$x, y, method = "sir", ndr = 1, keep.data = TRUE,
    maxitr = 40, verbose = FALSE
  )
  pred1 <- predict(fit, d$testx)
  pred2 <- predict(fit, d$testx)
  expect_equal(pred1$pred, pred2$pred)
})

# ── 6. control parameters ────────────────────────────────────────────────────

test_that("orthoDr_reg accepts custom control parameters", {
  set.seed(1123)
  d <- make_regression_data(n = 80, p = 4)
  y <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(80, sd = 0.2)

  fit <- orthoDr_reg(
    d$x, y, method = "sir", ndr = 1, maxitr = 30,
    control = list(rho = 1e-3, tau = 1e-2, ftol = 1e-4),
    verbose = FALSE
  )
  expect_equal(class(fit), c("orthoDr", "fit", "reg"))
})

# ── 7. ncore ──────────────────────────────────────────────────────────────────

test_that("orthoDr_reg runs with ncore > 0", {
  set.seed(1234)
  d <- make_regression_data(n = 80, p = 4)
  y <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(80, sd = 0.2)

  fit <- orthoDr_reg(
    d$x, y, method = "sir", ndr = 1,
    keep.data = TRUE, maxitr = 30, ncore = 2, verbose = FALSE
  )
  expect_equal(class(fit), c("orthoDr", "fit", "reg"))
  expect_true(all(is.finite(fit$B)))

  pred <- predict(fit, d$testx)
  expect_true(all(is.finite(pred$pred)))
})

# ── 8. verbose output ────────────────────────────────────────────────────────

test_that("verbose produces iteration output", {
  set.seed(1345)
  d <- make_regression_data(n = 60, p = 3)
  y <- 1 + d$x[, 1] + rnorm(60)

  expect_output(
    orthoDr_reg(d$x, y, method = "sir", ndr = 1, maxitr = 10, verbose = 1),
    "secs"
  )
})

# ── 9. Edge cases ────────────────────────────────────────────────────────────

test_that("orthoDr_reg works with minimum p = 2", {
  set.seed(1456)
  n <- 80; p <- 2
  x <- matrix(rnorm(n * p), n, p)
  y <- 1 + 0.5 * x[, 1] + rnorm(n, sd = 0.2)
  testx <- matrix(rnorm(10 * p), 10, p)

  fit <- orthoDr_reg(
    x, y, method = "sir", ndr = 1,
    keep.data = TRUE, maxitr = 30, verbose = FALSE
  )
  expect_equal(class(fit), c("orthoDr", "fit", "reg"))
  expect_equal(dim(fit$B), c(p, 1L))

  pred <- predict(fit, testx)
  expect_equal(length(pred$pred), 10)
})

test_that("orthoDr_reg handles n < p with B.initial", {
  # initB uses dr() which requires n > p, but we can bypass it with B.initial
  set.seed(1567)
  n <- 30
  p <- 50
  x <- matrix(rnorm(n * p), n, p)
  y <- 1 + x[, 1] - 0.5 * x[, 2] + rnorm(n, sd = 0.3)

  B0 <- matrix(0, p, 1)
  B0[1, 1] <- 1

  fit <- orthoDr_reg(
    x, y, method = "sir", ndr = 1,
    B.initial = B0, maxitr = 30, verbose = FALSE
  )
  expect_equal(class(fit), c("orthoDr", "fit", "reg"))
  expect_equal(dim(fit$B), c(p, 1L))
})

test_that("orthoDr_reg handles small n", {
  set.seed(1678)
  n <- 20
  p <- 3
  x <- matrix(rnorm(n * p), n, p)
  y <- 1 + x[, 1] + rnorm(n)

  fit <- orthoDr_reg(
    x, y, method = "sir", ndr = 1, maxitr = 30, verbose = FALSE
  )
  expect_equal(class(fit), c("orthoDr", "fit", "reg"))
  expect_equal(dim(fit$B), c(p, 1L))
})

test_that("orthoDr_reg rejects non-matrix x", {
  expect_error(
    orthoDr_reg(as.data.frame(matrix(1:9, 3, 3)), 1:3),
    "x must be a matrix"
  )
  # letters[1:9] is character, not a matrix, so is.matrix check fires first
  expect_error(
    orthoDr_reg(letters[1:9], 1:9),
    "x must be a matrix"
  )
})

test_that("orthoDr_reg rejects mismatched dimensions", {
  x <- matrix(rnorm(30), 10, 3)
  expect_error(
    orthoDr_reg(x, 1:9, method = "sir"),
    "Number of observations do not match"
  )
})

test_that("orthoDr_reg caps ndr to ncol(x)", {
  set.seed(1890)
  x <- matrix(rnorm(60), 20, 3)
  y <- 1 + x[, 1] + rnorm(20)

  # ndr=10 is larger than ncol(x)=3 — should cap to 3 with a warning
  expect_warning(
    fit <- orthoDr_reg(
      x, y, method = "sir", ndr = 10, maxitr = 30, verbose = FALSE
    ),
    "ndr > 4 is not recommended"
  )
  expect_equal(ncol(fit$B), 3L)
})

# ── 10. Direction recovery verification ──────────────────────────────────────

test_that("SIR recovers linear direction under linear model", {
  set.seed(2012)
  n <- 200
  p <- 4
  x <- matrix(rnorm(n * p), n, p)
  true_dir <- c(1, -0.5, 0.3, 0)  # 4th variable is noise
  true_dir <- true_dir / sqrt(sum(true_dir^2))
  y <- as.numeric(x %*% true_dir) + rnorm(n, sd = 0.2)

  fit <- orthoDr_reg(
    x, y, method = "sir", ndr = 1, keep.data = TRUE,
    maxitr = 80, verbose = FALSE
  )

  # Absolute correlation between estimated and true direction
  abs_cor <- abs(cor(fit$B[, 1], true_dir))
  expect_true(abs_cor > 0.9)
})

test_that("different methods can produce different estimates", {
  # SIR models the mean, PHD models the variance.
  # With a strong quadratic signal, the estimated directions should differ.
  set.seed(2123)
  n <- 100
  p <- 4
  x <- matrix(rnorm(n * p), n, p)
  # Strong quadratic signal with noise
  y <- (x[, 1]^2 - 0.6 * x[, 2]^2) * 3 + rnorm(n, sd = 0.15)

  fit_sir <- orthoDr_reg(
    x, y, method = "sir", ndr = 1, maxitr = 60, verbose = FALSE
  )
  fit_phd <- orthoDr_reg(
    x, y, method = "phd", ndr = 1, maxitr = 60, verbose = FALSE
  )
  # SIR targets mean (small for quadratic), PHD targets variance (large)
  # They should produce different objective values
  expect_true(abs(fit_sir$fn - fit_phd$fn) > 0.01)
  # Both should be valid
  expect_true(fit_sir$converge)
  expect_true(fit_phd$converge)
})

# ── 11. hMave ─────────────────────────────────────────────────────────────────

test_that("hMave returns expected structure", {
  set.seed(2345)
  n <- 100
  p <- 4
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  linpred <- 0.4 * x[, 1] - 0.2 * x[, 2]
  t_event <- rexp(n, rate = exp(linpred))
  t_censor <- rexp(n, rate = 0.8)
  y <- as.matrix(pmin(t_event, t_censor))
  censor <- as.matrix(as.numeric(t_event <= t_censor))

  fit <- hMave(x, y, censor, m0 = 1)

  expect_true(is.list(fit))
  expect_true(all(c("B", "cv") %in% names(fit)))
  expect_equal(dim(fit$B), c(p, 1L))
  expect_true(is.finite(fit$cv))
  expect_true(all(is.finite(fit$B)))
})

test_that("hMave works with custom B0", {
  set.seed(2456)
  n <- 80
  p <- 4
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  linpred <- 0.4 * x[, 1] - 0.2 * x[, 2]
  t_event <- rexp(n, rate = exp(linpred))
  t_censor <- rexp(n, rate = 1)
  y <- as.matrix(pmin(t_event, t_censor))
  censor <- as.matrix(as.numeric(t_event <= t_censor))

  B0 <- matrix(c(1, 0, 0, 0), 4, 1)
  fit <- hMave(x, y, censor, m0 = 1, B0 = B0)

  expect_equal(dim(fit$B), c(p, 1L))
  expect_true(all(is.finite(fit$B)))
})

test_that("hMave validates inputs", {
  n <- 60; p <- 3
  x <- matrix(rnorm(n * p), n, p)
  y <- as.matrix(rexp(n))
  censor <- as.matrix(rbinom(n, 1, 0.8))

  expect_error(hMave(as.data.frame(x), y, censor, 1), "x must be a matrix")
  expect_error(hMave(x, as.numeric(y), censor, 1), "y must be a matrix")
  expect_error(hMave(x, y, as.numeric(censor), 1), "censor must be a matrix")

  bad_censor <- as.matrix(c(rep(0, 30), rep(2, 30)))
  expect_error(hMave(x, y, bad_censor, 1), "censor can only have two levels 0/1")

  # Bad B0: wrong nrow (4 rows for 3-column x)
  bad_B0_nrow <- matrix(1:4, 4, 1)
  expect_error(hMave(x, y, censor, 1, B0 = bad_B0_nrow), "Dimension of B0 is not correct")
})

# ── 12. pSAVE ─────────────────────────────────────────────────────────────────

test_that("pSAVE returns expected structure", {
  set.seed(2567)
  n <- 200; p <- 4
  x <- matrix(rnorm(n * p), n, p)
  a <- runif(n, 0, 2)
  # Simpler signal: linear in first two columns
  r <- 1 + 0.5 * x[, 1] - 0.3 * x[, 2] + 0.2 * a + rnorm(n, sd = 0.5)

  B <- pSAVE(x, a, r, ndr = 2)

  expect_true(is.matrix(B))
  expect_equal(dim(B), c(p, 2L))
  expect_true(all(is.finite(B)))
})

test_that("pSAVE validates inputs", {
  x <- matrix(rnorm(30), 10, 3)
  expect_error(pSAVE(as.data.frame(x), 1:10, 1:10), "X must be a matrix")
  expect_error(pSAVE(x, 1:9, 1:10), "Number of observations do not match")
})

# ── 13. Class chain and print method ─────────────────────────────────────────

test_that("print.orthoDr works on fitted regression", {
  set.seed(2678)
  d <- make_regression_data(n = 60, p = 3)
  y <- 1 + d$x[, 1] + rnorm(60)

  fit <- orthoDr_reg(
    d$x, y, method = "sir", ndr = 1, maxitr = 30, verbose = FALSE
  )
  expect_output(print(fit), "Subspace for reg model using sir")
})

test_that("print.orthoDr works on regression prediction", {
  set.seed(2789)
  d <- make_regression_data(n = 60, p = 3)
  y <- 1 + d$x[, 1] + rnorm(60)

  fit <- orthoDr_reg(
    d$x, y, method = "sir", ndr = 1, keep.data = TRUE,
    maxitr = 30, verbose = FALSE
  )
  pred <- predict(fit, d$testx)
  expect_output(print(pred), "Prediction for orthoDr regression")
})

# ── 14. Reproducibility ──────────────────────────────────────────────────────

test_that("orthoDr_reg is reproducible with set.seed", {
  set.seed(2890)
  x <- matrix(rnorm(80), 40, 2)
  y <- 1 + x[, 1] + rnorm(40)

  set.seed(42)
  fit1 <- orthoDr_reg(
    x, y, method = "sir", ndr = 1, maxitr = 50, verbose = FALSE
  )
  set.seed(42)
  fit2 <- orthoDr_reg(
    x, y, method = "sir", ndr = 1, maxitr = 50, verbose = FALSE
  )
  expect_equal(fit1$B, fit2$B)
  expect_equal(fit1$fn, fit2$fn)
})

# ── 15. Numerical stability ──────────────────────────────────────────────────

test_that("orthoDr_reg returns orthonormal columns", {
  set.seed(3001)
  d <- make_regression_data(n = 100, p = 5)
  y <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(100, sd = 0.3)

  for (ndr_val in c(1, 2)) {
    fit <- orthoDr_reg(
      d$x, y, method = "sir", ndr = ndr_val,
      maxitr = 50, verbose = FALSE
    )
    BtB <- crossprod(fit$B)
    # Diagonal should be 1
    expect_equal(diag(as.matrix(BtB)), rep(1, ndr_val), tolerance = 1e-5)
    # Off-diagonal should be near 0
    off <- BtB - diag(ndr_val)
    expect_equal(max(abs(off)), 0, tolerance = 1e-5)
  }
})

test_that("prediction handles new data with different nrow", {
  set.seed(3112)
  d <- make_regression_data(n = 100, p = 4)
  y <- 1 + d$x[, 1] - 0.5 * d$x[, 2] + rnorm(100, sd = 0.3)

  fit <- orthoDr_reg(
    d$x, y, method = "sir", ndr = 1, keep.data = TRUE,
    maxitr = 50, verbose = FALSE
  )

  # Predict for 1, 5, and 50 new observations
  for (ntest in c(1, 5, 50)) {
    tx <- matrix(rnorm(ntest * 4), ntest, 4)
    pred <- predict(fit, tx)
    expect_equal(length(pred$pred), ntest)
    expect_true(all(is.finite(pred$pred)))
  }
})
