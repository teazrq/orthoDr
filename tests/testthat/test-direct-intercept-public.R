intercept_public_fixture <- function() {
  set.seed(48213)
  x <- matrix(runif(60, -1, 1), 20, 3)
  a <- seq(.05, 1.95, length.out = 20)
  r <- 3 + .2 * x[, 1] - (a - 1 - .3 * x[, 2])^2
  list(x = x, a = a, r = r, ndr = 1L,
       B.initial = matrix(c(1, 2, -1) / sqrt(6), 3, 1),
       method = "direct", lambda = .1, lambda.reward = .1,
       K = 8L, bw = .7, maxitr = 0L, keep.data = TRUE, ncore = 1L)
}

intercept_public_fit <- function(...) {
  do.call(orthoDr_pdose, modifyList(intercept_public_fixture(), list(...)))
}

# Independent cross-kernel construction uses the stated Gaussian equation.
intercept_public_kernel <- function(fit, newx) {
  train <- fit$x %*% fit$B
  test <- newx %*% fit$B
  scale <- apply(train, 2, sd)
  exponent <- matrix(0, nrow(newx), nrow(train))
  for (j in seq_len(ncol(train)))
    exponent <- exponent - outer(test[, j], train[, j], "-")^2 /
      (2 * fit$bw^2 * scale[j]^2)
  exp(exponent)
}

test_that("the intercept API validates only its intended direct mode", {
  for (bad in list(NA, 1, "TRUE", logical(), c(TRUE, FALSE)))
    expect_error(intercept_public_fit(intercept = bad), "intercept must")
  for (bad in list(c(0, 0), c(2, 0), 1, c(0, Inf), c(NA, 2),
                   c(0, 1 + 1i), matrix(c(0, 2), 1), "0,2"))
    expect_error(intercept_public_fit(intercept = TRUE, intercept.bounds = bad),
                 "intercept.bounds")
  expect_error(intercept_public_fit(intercept.bounds = c(0, 2)), "requires intercept")
  expect_error(intercept_public_fit(intercept = TRUE, lambda.reward = 0),
               "lambda.reward must be positive")
  expect_no_error(intercept_public_fit(lambda.reward = 0))

  pseudo <- intercept_public_fixture()
  pseudo$method <- "pseudo_direct"
  pseudo$lambda.reward <- NULL
  expect_error(do.call(orthoDr_pdose, c(pseudo, list(intercept = FALSE))),
               "only used with method")
  expect_error(do.call(orthoDr_pdose, c(pseudo, list(intercept.bounds = NULL))),
               "only used with method")
})

test_that("new fits retain centered prediction state and a training mean intercept", {
  fit <- intercept_public_fit(intercept = TRUE)
  expect_identical(fit$objective.version, "direct_joint_intercept_v2")
  expect_equal(as.vector(fit$intercept.bounds), range(fit$a))
  train_kernel <- intercept_public_kernel(fit, fit$x)
  expect_equal(as.vector(fit$kernel.mean), colMeans(train_kernel), tolerance = 1e-12)
  expect_equal(mean(predict(fit, fit$x)$pred), fit$intercept, tolerance = 1e-12)
  newx <- rbind(fit$x[1:3, , drop = FALSE], c(1e4, -1e4, 1e4))
  expected <- fit$intercept + drop(sweep(intercept_public_kernel(fit, newx),
                                        2, as.vector(fit$kernel.mean), "-") %*% fit$W)
  expect_equal(as.vector(predict(fit, newx)$pred), expected, tolerance = 1e-12)
  expect_equal(tail(expected, 1), fit$intercept - sum(fit$kernel.mean * fit$W),
               tolerance = 1e-12)
  expect_named(fit$gradient.norm, c("B", "W", "intercept", "policy"))
  expect_equal(unname(fit$gradient.norm["policy"]),
               sqrt(sum(fit$gradient.norm[c("W", "intercept")]^2)), tolerance = 1e-12)

  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  saveRDS(fit, saved)
  expect_identical(predict(readRDS(saved), newx), predict(fit, newx))
  original_intercept <- fit$intercept
  fit$intercept <- 0
  expect_equal(as.vector(predict(fit, newx)$pred), expected - original_intercept,
               tolerance = 1e-12)
})

test_that("intercept bounds constrain the mean without clipping predictions", {
  fit <- intercept_public_fit(intercept = TRUE, intercept.bounds = c(-3, -2))
  expect_equal(fit$intercept, -2)
  expect_equal(mean(predict(fit, fit$x)$pred), -2, tolerance = 1e-12)
  expect_error(intercept_public_fit(intercept = TRUE, intercept.bounds = c(-3, -2),
                                    lambda = c(.1, .2)), "scalar lambda")

  # A stored finite policy can leave the mean-dose interval for individual rows.
  fit$W <- as.vector(fit$W) + 20 * seq_along(fit$W)
  predictions <- as.vector(predict(fit, fit$x)$pred)
  expect_equal(mean(predictions), fit$intercept, tolerance = 1e-9)
  expect_true(any(predictions < -3 | predictions > -2))
})

test_that("constant grid targets have zero coefficients including at a bound", {
  fit <- intercept_public_fit(intercept = TRUE, K = 1, lambda = c(.3, .1))
  expect_equal(fit$intercept, min(fit$a))
  expect_equal(as.vector(fit$W), rep(0, nrow(fit$x)))
  expect_equal(fit$lambda.init, .3)
  expect_equal(as.vector(predict(fit, fit$x)$pred), rep(min(fit$a), nrow(fit$x)))
})

test_that("legacy predictions and data-retention behavior are preserved", {
  legacy <- intercept_public_fit()
  expect_identical(legacy, intercept_public_fit(intercept = FALSE))
  expect_null(legacy$intercept)
  expect_null(legacy$kernel.mean)
  expected <- predict(legacy, legacy$x)
  legacy$objective.version <- NULL
  expect_identical(predict(legacy, legacy$x), expected)
  expect_equal(as.vector(expected$pred),
               drop(intercept_public_kernel(legacy, legacy$x) %*% legacy$W), tolerance = 1e-12)

  hidden <- intercept_public_fit(intercept = TRUE, keep.data = FALSE)
  expect_length(hidden$kernel.mean, nrow(intercept_public_fixture()$x))
  expect_true(is.finite(hidden$intercept))
  expect_error(predict(hidden, intercept_public_fixture()$x), "keep.data = TRUE")
})

test_that("incomplete new prediction state is not interpreted as a legacy fit", {
  fit <- intercept_public_fit(intercept = TRUE)
  for (field in c("intercept", "kernel.mean")) {
    broken <- fit
    broken[[field]] <- NULL
    expect_error(predict(broken, fit$x), "missing valid intercept or kernel.mean")
    broken[[field]] <- NA_real_
    expect_error(predict(broken, fit$x), "missing valid intercept or kernel.mean")
  }
  broken <- fit
  broken$objective.version <- NULL
  expect_error(predict(broken, fit$x), "supported objective.version")
})
