test_that("orthoDr_surv supports dm, dn, and forward", {
  set.seed(456)
  n <- 100
  p <- 4
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)

  linpred <- 0.5 * x[, 1] - 0.25 * x[, 2]
  t_event <- rexp(n, rate = exp(linpred))
  t_censor <- rexp(n, rate = 0.8)
  y <- pmin(t_event, t_censor)
  censor <- as.numeric(t_event <= t_censor)
  testx <- matrix(rnorm(12 * p), nrow = 12, ncol = p)

  methods <- c("dm", "dn", "forward")

  for (method in methods) {
    fit <- orthoDr_surv(
      x,
      y,
      censor,
      method = method,
      ndr = 1,
      keep.data = TRUE,
      maxitr = 40,
      verbose = FALSE
    )

    expect_equal(class(fit), c("orthoDr", "fit", "surv"))
    expect_equal(dim(fit$B), c(p, 1))
    expect_true(all(is.finite(fit$B)))

    if (method %in% c("dm", "dn")) {
      pred <- predict(fit, testx)

      expect_equal(class(pred), c("orthoDr", "predict", "surv"))
      expect_true(is.matrix(pred$surv))
      expect_equal(ncol(pred$surv), nrow(testx))
      expect_true(length(pred$timepoints) >= 1)
      expect_true(all(is.finite(pred$surv)))
      expect_true(all(pred$surv >= 0 & pred$surv <= 1))
    }
  }
})

test_that("CP_SIR returns expected structure", {
  set.seed(567)
  n <- 90
  p <- 4
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  linpred <- 0.4 * x[, 1] - 0.3 * x[, 2]
  t_event <- rexp(n, rate = exp(linpred))
  t_censor <- rexp(n, rate = 0.9)
  y <- pmin(t_event, t_censor)
  censor <- as.numeric(t_event <= t_censor)

  fit_cpsir <- CP_SIR(x, y, censor)

  expect_true(is.list(fit_cpsir))
  expect_true(all(c("values", "vectors") %in% names(fit_cpsir)))
  expect_equal(length(fit_cpsir$values), p)
  expect_equal(dim(fit_cpsir$vectors), c(p, p))
  expect_true(all(is.finite(fit_cpsir$values)))
  expect_true(all(is.finite(fit_cpsir$vectors)))
})

test_that("survival prediction requires keep.data", {
  set.seed(654)
  n <- 80
  p <- 3
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  t_event <- rexp(n, rate = exp(0.3 * x[, 1]))
  t_censor <- rexp(n, rate = 1)
  y <- pmin(t_event, t_censor)
  censor <- as.numeric(t_event <= t_censor)

  fit <- orthoDr_surv(
    x,
    y,
    censor,
    method = "dm",
    ndr = 1,
    keep.data = FALSE,
    maxitr = 30,
    verbose = FALSE
  )

  testx <- matrix(rnorm(6 * p), nrow = 6, ncol = p)
  expect_error(
    predict(fit, testx),
    "Need the original data for prediction"
  )
})