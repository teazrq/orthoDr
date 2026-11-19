reg_coordinate_objective <- function(x, Y, B, bw, method) {
  weights <- function(z) {
    kernel <- exp(-as.matrix(dist(z))^2)
    kernel / rowSums(kernel)
  }
  centered <- scale(x, scale = FALSE)
  projected <- centered %*% B
  wx <- weights(scale(projected, center = FALSE, scale = apply(projected, 2, sd) * bw * sqrt(2)))
  wy <- weights(Y)
  if (method == "sir") {
    ey <- wy %*% centered
    moment <- crossprod(ey - wx %*% ey, centered - wx %*% centered) / nrow(x)
  } else {
    local_x <- wx %*% centered
    local_y <- wy %*% centered
    moment <- matrix(0, ncol(x), ncol(x))
    for (i in seq_len(nrow(x))) {
      centered_x <- sweep(centered, 2, local_x[i, ], "-")
      centered_y <- sweep(centered, 2, local_y[i, ], "-")
      cx <- crossprod(centered_x, centered_x * wx[i, ])
      cy <- crossprod(centered_y, centered_y * wy[i, ])
      moment <- moment + (diag(ncol(x)) - cy) %*%
        (tcrossprod(centered_y[i, ], centered_x[i, ]) - cx)
    }
    moment <- moment / nrow(x)
  }
  sum(moment^2)
}

test_that("SIR and SAVE match fixed-basis moments and tangent derivatives", {
  set.seed(139)
  x <- matrix(rnorm(128), 32, 4)
  y <- x[, 1] + .4 * x[, 2]^2 + rnorm(nrow(x), sd = .2)
  Y <- scale(y) * length(y)^(1 / 6) / sqrt(2)
  bw <- .7
  base <- diag(4)[, 1, drop = FALSE]
  tangent <- diag(4)[, 2, drop = FALSE]
  h <- 1e-5
  bases <- list(base, cos(h) * base + sin(h) * tangent, cos(h) * base - sin(h) * tangent)
  for (method in c("sir", "save")) {
    ref <- vapply(bases, function(B) reg_coordinate_objective(x, Y, B, bw, method), 0.0)
    transformed <- if (method == "sir") .01 * x else 1e5 + x
    conversion <- if (method == "sir") .01^4 else 1
    native <- vapply(bases, function(B) {
      orthoDr:::reg_init(method, B, transformed, Y, bw, 1L) / conversion
    }, 0.0)
    expect_equal(native, ref, tolerance = 2e-10)
    expect_equal((native[2] - native[3]) / (2 * h),
                 (ref[2] - ref[3]) / (2 * h), tolerance = 2e-7)
    fit <- orthoDr_reg(transformed, y, method, 1, base, bw,
                       control = list(gtol = 1e5), ncore = 1)
    expect_equal(fit$itr, 0)
    expect_equal(fit$fn / conversion, ref[1], tolerance = 2e-10)
  }
})

test_that("SAVE anchors both local covariance calculations", {
  groups_x <- c(0, 1, 1, 0, 1, 0, 1, 1, 3, 3, 2, 2, 2, 2, 3, 2)
  groups_y <- c(1, 0, 0, 0, 0, 0, 0, 1, 2, 3, 3, 3, 2, 3, 3, 2)
  local_x <- cbind(groups_x, rep(c(rep(0, 7), 1), 2))
  x <- local_x + cbind(0, rep(c(-2^50, 2^50), each = 8))
  Y <- matrix(100 * groups_y)
  B <- matrix(c(1, 0))
  bw <- .005
  # Matching Gaussian labels have weight 1. Each group has a common offset,
  # so its covariance can be independently calculated in small coordinates.
  moment <- matrix(0, 2, 2)
  for (i in seq_len(nrow(x))) {
    x_group <- local_x[groups_x == groups_x[i], , drop = FALSE]
    y_group <- local_x[groups_y == groups_y[i], , drop = FALSE]
    rx <- local_x[i, ] - colMeans(x_group)
    ry <- local_x[i, ] - colMeans(y_group)
    cx <- crossprod(scale(x_group, scale = FALSE)) / nrow(x_group)
    cy <- crossprod(scale(y_group, scale = FALSE)) / nrow(y_group)
    moment <- moment + (diag(2) - cy) %*% (tcrossprod(ry, rx) - cx)
  }
  reference <- sum((moment / nrow(x))^2)
  expect_equal(orthoDr:::save_init(B, x, Y, bw, 1L), reference, tolerance = 1e-14)
  fit <- orthoDr:::save_solver(B, x, Y, bw, 1e-4, .2, .85, 1e-3,
                              1e-6, 1e-6, 1e-6, 1e100, 5L, 0L, 1L)
  expect_equal(fit$itr, 0)
  expect_equal(fit$fn, reference, tolerance = 1e-14)
})

test_that("small predictor units do not halt regression at its initial directions", {
  set.seed(117)
  x <- matrix(round(32 * rnorm(224)) / 32, 56, 4)
  y <- x[, 1]^2 + .5 * x[, 2] + rnorm(nrow(x), sd = .15)
  factor <- 2^-10
  for (method in c("phd", "local", "seff")) {
    power <- if (method == "phd") 4L else 2L
    for (start in list(NULL, matrix(c(.3, 1, .5, -.25)))) {
      run <- function(x) orthoDr_reg(x, y, method, 1L, start, .65,
                                     keep.data = TRUE, maxitr = 12L, ncore = 1L)
      base <- run(x)
      scaled <- run(factor * x)
      expect_gt(scaled$itr, 0)
      expect_equal(scaled$itr, base$itr)
      expect_equal(tcrossprod(scaled$B), tcrossprod(base$B), tolerance = 1e-10)
      expect_equal(scaled$fn / factor^power, base$fn, tolerance = 1e-10)
      query <- x[2:6, , drop = FALSE] + .125
      expect_equal(predict(scaled, factor * query), predict(base, query), tolerance = 1e-10)
    }
  }
})

test_that("PHD uses centered outer products for its reported objective", {
  set.seed(117)
  x <- matrix(round(32 * rnorm(224)) / 32, 56, 4)
  y <- x[, 1]^2 + .5 * x[, 2] + rnorm(nrow(x), sd = .15)
  shifted <- sweep(x, 2, c(10, 3, -2, 5), "+")
  for (start in list(NULL, matrix(c(.3, 1, .5, -.25)))) {
    run <- function(x) orthoDr_reg(x, y, "phd", 1L, start, .65,
                                   keep.data = TRUE, maxitr = 12L, ncore = 1L)
    base <- run(x)
    fit <- run(shifted)
    expect_equal(tcrossprod(fit$B), tcrossprod(base$B), tolerance = 1e-10)
    Y <- scale(y) * length(y)^(1 / 6) / sqrt(2)
    X <- scale(x, scale = FALSE)
    z <- X %*% fit$B
    W <- exp(-(outer(z[, 1], z[, 1], "-") / (sd(z) * .65 * sqrt(2)))^2)
    W <- W / rowSums(W)
    yy <- Y - W %*% Y
    moment <- matrix(0, ncol(x), ncol(x))
    for (j in seq_len(ncol(x))) {
      xx <- X * X[, j]
      moment[, j] <- colMeans((xx - W %*% xx) * as.vector(yy))
    }
    expect_equal(fit$fn, sum(moment^2), tolerance = 1e-10)
    expect_equal(predict(fit, shifted[2:6, , drop = FALSE] + .125),
                 predict(base, x[2:6, , drop = FALSE] + .125), tolerance = 1e-10)
  }
})

test_that("SAVE descends at small scales and reports its supplied-unit objective", {
  set.seed(139)
  x <- matrix(rnorm(128), 32, 4)
  y <- x[, 1] + .4 * x[, 2]^2 + rnorm(nrow(x), sd = .2)
  Y <- scale(y) * length(y)^(1 / 6) / sqrt(2)
  for (factor in c(.001, 16)) {
    for (start in list(NULL, matrix(c(.3, 1, .5, -.25)))) {
      run <- function(gtol) orthoDr_reg(factor * x, y, "save", 1L, start, .7,
                                        maxitr = 12L, control = list(gtol = gtol), ncore = 1L)
      initial <- run(1e100)
      fit <- run(1e-6)
      expect_equal(initial$itr, 0)
      expect_gt(fit$itr, 0)
      expect_lt(fit$fn, initial$fn)
      reference <- reg_coordinate_objective(factor * x, Y, fit$B, .7, "save")
      expect_equal(fit$fn / reference, 1, tolerance = 1e-10)
    }
  }
})

test_that("SIR can report a finite objective after an extreme change of units", {
  x <- matrix(c(-1, 1))
  y <- c(0, 1)
  query <- matrix(c(0, .5))
  base <- orthoDr_reg(x, y, keep.data = TRUE, ncore = 1)
  for (power in c(256, -260)) {
    factor <- 2^power
    fit <- orthoDr_reg(factor * x, y, keep.data = TRUE, ncore = 1)
    expected <- (base$fn * 2^(2 * power)) * 2^(2 * power)
    expect_true(is.finite(fit$fn) && fit$fn > 0)
    expect_equal(fit$fn / expected, 1, tolerance = 1e-11)
    expect_identical(fit$B, base$B)
    expect_equal(predict(fit, factor * query), predict(base, query), tolerance = 1e-14)
  }
  factor <- .Machine$double.xmax
  fit <- orthoDr_reg(factor * x, y, keep.data = TRUE, ncore = 1)
  expect_true(fit$converge)
  expect_identical(fit$fn, Inf)
  expect_equal(predict(fit, factor * query), predict(base, query), tolerance = 1e-14)
})

test_that("narrow SIR and SAVE kernels retain close observations after a distant row", {
  t <- 0:63
  x <- cbind(c(2^54, t %% 4), c(0, round(16 * sin(1.21 * t)) / 16))
  bw <- .35 / sd(x[, 1])
  for (method in c("sir", "save")) {
    y <- c(if (method == "sir") 10 else 1e6, sin(2.1 * t))
    reference <- if (method == "sir") 1.89117288877e-10 else 9.2609525559169e-6
    for (order in list(1:65, c(2:65, 1))) {
      fit <- orthoDr_reg(x[order, ], y[order], method, 1, matrix(c(1, 0)), bw,
                         control = list(gtol = 1e100), ncore = 1)
      expect_equal(fit$itr, 0)
      expect_equal(fit$fn / reference, 1, tolerance = 1e-8)
    }
  }
})

test_that("SIR retains narrow neighborhoods in both separated clusters", {
  t <- 0:63
  local <- rbind(cbind(.25 * (t %% 4), round(16 * sin(1.21 * t)) / 16),
                 cbind(.25 * (t %% 4), round(16 * sin(1.21 * t)) / 16))
  groups <- rep(c(-1, 1), each = length(t))
  offset <- 2^50
  x <- local
  x[, 1] <- x[, 1] + offset * groups
  y <- 2 * groups + rep(sin(2.1 * t), 2)
  bw <- .35 * .25 / sd(x[, 1])
  Y <- as.vector(scale(y)) * nrow(x)^(1 / 6) / sqrt(2)
  wx <- exp(-(outer(local[, 1], local[, 1], "-") / (.35 * .25 * sqrt(2)))^2)
  wx <- wx * outer(groups, groups, "==")
  wx <- wx / rowSums(wx)
  wy <- exp(-outer(Y, Y, "-")^2)
  wy <- wy / rowSums(wy)
  # Calculate the responses' conditional mean for each cluster offset
  # independently of the within-cluster coordinates.
  ey <- wy %*% local
  residual_ey <- ey - wx %*% ey
  ey_group <- wy %*% groups
  residual_ey[, 1] <- residual_ey[, 1] + offset * (ey_group - wx %*% ey_group)
  moment <- crossprod(residual_ey, local - wx %*% local) / nrow(x)
  reference <- sum(moment^2)
  for (order in list(1:128, c(65:128, 1:64))) {
    fit <- orthoDr_reg(x[order, ], y[order], ndr = 1, B.initial = matrix(c(1, 0)),
                       bw = bw, control = list(gtol = 1e100), ncore = 1)
    expect_equal(fit$itr, 0)
    expect_equal(fit$fn / reference, 1, tolerance = 1e-10)
  }
})

test_that("SIR and SAVE retain a contrast between large nearly equal columns", {
  x <- rbind(c(2^54, 2^54), cbind(rep(c(0, 1, 0, 1), 3), rep(c(0, 0, 1, 1), 3)))
  y <- sin(1.7 * (0:12))
  B <- matrix(c(1, -1)) / sqrt(2)
  for (method in c("sir", "save")) {
    reference <- if (method == "sir") 2.4256406167671955e61 else 2.131848331387999e126
    for (order in list(1:13, c(2:13, 1))) {
      fit <- orthoDr_reg(x[order, ], y[order], method, 1, B, bw = .6,
                         control = list(gtol = 1e150), ncore = 1)
      expect_equal(fit$itr, 0)
      expect_equal(fit$fn / reference, 1, tolerance = 1e-10)
    }
  }
})

test_that("supplied SIR and SAVE starts can bypass sliced fitting when p exceeds n", {
  x <- matrix(0:19, 4, 5, byrow = TRUE)
  y <- c(0, 1, 3, 2)
  for (method in c("sir", "save")) {
    fit <- orthoDr_reg(x, y, method, 1, diag(5)[, 1, drop = FALSE],
                       maxitr = 5, keep.data = TRUE, ncore = 1)
    expect_true(is.finite(fit$fn))
    expect_true(all(is.finite(predict(fit, x)$pred)))
  }
})
