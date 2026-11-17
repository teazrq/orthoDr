# These references use the mathematical pair derivatives and explicit smoother
# matrices, independently of the native vectorized implementation.
pseudo_native_fixture <- function(p = 6L, d = 2L) {
  set.seed(781L + p + d)
  n <- 13L
  x <- matrix(rnorm(n * p), n, p)
  x <- sweep(x, 2, c(1, 4, .5, 2, 3, .8)[seq_len(p)], "*")
  x[, 2] <- x[, 2] + .6 * x[, 1]
  B <- qr.Q(qr(matrix(sin(seq_len(p * d)), p, d)))
  a <- seq(.1, 3.8, length.out = n)
  r <- 2 * exp(-(a - 2 - .5 * tanh(x[, 1]))^2) + .1 * x[, 3]
  list(x = x, B = B, a = a, r = r, bw = .73,
       grid = seq(min(a), max(a), length.out = 9L), lambda = c(.03, .2, 1))
}

pseudo_native_reference <- function(z, B = z$B, intercept = FALSE,
                                    dose.factor = 1 / sqrt(2)) {
  n <- nrow(z$x)
  q <- (z$r - mean(z$r)) / sd(z$r)
  u <- z$x %*% B
  s2 <- apply(u, 2, var)
  logk <- matrix(0, n, n)
  for (k in seq_len(ncol(B)))
    logk <- logk - outer(u[, k], u[, k], "-")^2 / (2 * z$bw^2 * s2[k])
  loga <- -outer(z$a, z$a, "-")^2 / (2 * z$bw^2 * var(z$a))
  stable <- function(logw) {
    w <- exp(logw - max(logw))
    w / sum(w)
  }
  probability <- apply(logk + loga, 2, stable)
  fitted <- drop(crossprod(probability, q))
  objective <- mean((q - fitted)^2)
  sensitivity <- sweep(probability * sweep(matrix(q, n, n), 2, fitted, "-"),
                        2, 2 * (fitted - q) / n, "*")
  gradient <- matrix(0, nrow(B), ncol(B))
  covariance <- cov(z$x)
  for (k in seq_len(ncol(B))) for (i in seq_len(n)) for (j in seq_len(n)) {
    delta <- u[j, k] - u[i, k]
    pair <- -delta * (z$x[j, ] - z$x[i, ]) / (z$bw^2 * s2[k]) +
      delta^2 * drop(covariance %*% B[, k]) / (z$bw^2 * s2[k]^2)
    gradient[, k] <- gradient[, k] + sensitivity[j, i] * pair
  }
  scores <- matrix(0, n, length(z$grid))
  for (i in seq_len(n)) for (j in seq_along(z$grid)) {
    logw <- logk[, i] - (z$a - z$grid[j])^2 /
      (2 * z$bw^2 * var(z$a) * dose.factor^2)
    scores[i, j] <- sum(stable(logw) * q)
  }
  targets <- z$grid[max.col(scores, ties.method = "first")]
  K <- exp(logk)
  ridges <- lapply(z$lambda, function(lambda) {
    Cinv <- solve(K + diag(lambda, n))
    if (intercept) {
      v <- drop(Cinv %*% rep(1, n)); t <- sum(v)
      alpha <- sum(v * targets) / t
      W <- drop(Cinv %*% targets) - v * alpha
      H <- diag(n) - lambda * Cinv + lambda * tcrossprod(v) / t
    } else {
      alpha <- 0
      W <- drop(Cinv %*% targets)
      H <- diag(n) - lambda * Cinv
    }
    gcv <- n * sum((targets - H %*% targets)^2) / (n - sum(diag(H)))^2
    list(W = W, alpha = alpha, GCV = gcv)
  })
  gcv <- vapply(ridges, `[[`, numeric(1), "GCV")
  selected <- which.min(gcv)
  list(objective = objective, gradient = gradient,
       residual = gradient - B %*% (t(gradient) %*% B),
       fitted = fitted, kernel = K, scores = scores, targets = targets,
       GCV = gcv, W = ridges[[selected]]$W,
       alpha = ridges[[selected]]$alpha, lambda = z$lambda[selected])
}

pseudo_native_eval <- function(z, B = z$B, intercept = FALSE, step = 1e-4,
                               factor = 1 / sqrt(2)) {
  orthoDr:::pdose_semi_evaluate(B, z$x, z$r, z$a, z$bw,
    z$grid, z$lambda, intercept, step, factor)
}

pseudo_native_audit <- function(z, maxitr = 0L, intercept = FALSE,
                                control = list(), budget = 30L,
                                fail.gradient = -1L, fail.trial = -1L) {
  ctl <- modifyList(list(rho = 1e-4, eta = .2, gamma = .85, tau = 1,
    btol = 1e-10, ftol = 1e-10, gtol = 1e-6), control)
  orthoDr:::pdose_semi_audit_solver(z$B, z$x, z$r, z$a, z$grid,
    z$lambda, z$bw, ctl$rho, ctl$eta, ctl$gamma, ctl$tau, 1e-6,
    ctl$btol, ctl$ftol, ctl$gtol, maxitr, 0L, 1L, intercept,
    TRUE, budget, fail.gradient, fail.trial)
}

pseudo_native_events <- function(fit, event) {
  Filter(function(item) identical(item$event, event), fit$trace)
}

test_that("normalized MSE and projected-variance derivatives match independent equations", {
  for (p in c(3L, 6L)) for (d in 1:2) {
    z <- pseudo_native_fixture(p, d)
    B <- z$B + matrix(seq_along(z$B), nrow(z$B)) / 200
    actual <- pseudo_native_eval(z, B)
    expected <- pseudo_native_reference(z, B)
    expect_equal(actual$objective.scaled, expected$objective, tolerance = 1e-12)
    expect_equal(actual$fn, sd(z$r) * sqrt(expected$objective), tolerance = 1e-12)
    expect_equal(actual$gradient, expected$gradient, tolerance = 1e-11)
    expect_equal(actual$residual, expected$residual, tolerance = 1e-11)
    eps <- 2e-6
    finite.difference <- matrix(0, nrow(B), ncol(B))
    for (k in seq_along(B)) {
      delta <- B * 0; delta[k] <- eps
      plus <- pseudo_native_eval(z, B + delta)$objective.scaled
      minus <- pseudo_native_eval(z, B - delta)$objective.scaled
      finite.difference[k] <- (plus - minus) / (2 * eps)
    }
    expect_equal(actual$gradient, finite.difference, tolerance = 1e-7)
    tangent <- matrix(cos(seq_along(B)), nrow(B))
    tangent <- tangent - z$B %*% ((crossprod(z$B, tangent) + crossprod(tangent, z$B)) / 2)
    for (h in c(1e-4, 1e-5, 1e-6)) {
      derivative <- (pseudo_native_eval(z, z$B + h * tangent)$objective.scaled -
        pseudo_native_eval(z, z$B - h * tangent)$objective.scaled) / (2 * h)
      native <- pseudo_native_eval(z)
      expect_lt(abs(sum(native$gradient * tangent) - derivative), 1e-7)
    }
  }
})

test_that("both Cayley branches use their actual curve slopes", {
  for (p in c(3L, 6L)) {
    z <- pseudo_native_fixture(p, 2L)
    e <- pseudo_native_eval(z, step = 1e-6)
    expected <- sum(e$gradient * e$residual)
    expect_equal(e$cayley.dense$slope, expected, tolerance = 1e-12)
    expect_equal(e$cayley.lowrank$slope, expected, tolerance = 1e-12)
    expect_true(e$cayley.dense$valid)
    expect_true(e$cayley.lowrank$valid)
    expect_equal(e$cayley.dense$B, e$cayley.lowrank$B, tolerance = 1e-12)
    expect_lt(norm(crossprod(e$cayley.dense$B) - diag(2), "F"), 1e-12)
    for (branch in c("cayley.dense", "cayley.lowrank")) {
      plus <- pseudo_native_eval(z, e[[branch]]$B)$objective.scaled
      expect_equal((plus - e$objective.scaled) / 1e-6,
                   -e[[branch]]$slope, tolerance = 1e-6)
    }
  }
})

test_that("reward affine changes preserve the objective and exact constants are rejected", {
  z <- pseudo_native_fixture()
  reference <- pseudo_native_eval(z)
  for (s in c(1e-8, 1, 1e8)) {
    transformed <- z; transformed$r <- s * (z$r + 3)
    actual <- pseudo_native_eval(transformed)
    expect_equal(actual$objective.scaled, reference$objective.scaled, tolerance = 1e-12)
    expect_equal(actual$gradient, reference$gradient, tolerance = 1e-11)
    expect_equal(actual$grid.scores, reference$grid.scores, tolerance = 1e-11)
    expect_equal(actual$fn, s * reference$fn, tolerance = 1e-12)
  }
  for (constant in c(.1, 0, -3)) {
    flat <- z; flat$r[] <- constant
    expect_error(pseudo_native_eval(flat), "Constant rewards.*identify")
    expect_error(pseudo_native_audit(flat), "Constant rewards.*identify")
  }
  z$bw <- .001
  exact <- pseudo_native_audit(z)
  expect_equal(exact$objective.scaled, 0, tolerance = 0)
  expect_equal(exact$gradient.norm, 0, tolerance = 0)
  expect_true(exact$converge)
  expect_identical(exact$stop.reason, "gradient_tolerance")
})

test_that("raw log-weight dose extraction is finite at narrow bandwidths", {
  for (bw in c(.05, .01, .001)) {
    z <- pseudo_native_fixture(); z$bw <- bw
    actual <- pseudo_native_eval(z)
    expected <- pseudo_native_reference(z)
    expect_true(all(is.finite(actual$grid.scores)))
    expect_equal(actual$grid.scores, expected$scores, tolerance = 1e-10)
    expect_equal(drop(actual$grid.targets), expected$targets, tolerance = 0)
    expect_equal(actual$dose.bandwidth.scale, sd(z$a) * bw / sqrt(2), tolerance = 1e-14)
  }
  z <- pseudo_native_fixture(); z$grid <- c(.1, .1)
  actual <- pseudo_native_eval(z)
  expect_equal(actual$grid.scores[, 1], actual$grid.scores[, 2], tolerance = 0)
  expect_equal(drop(actual$grid.targets), rep(.1, nrow(z$x)), tolerance = 0)
  aligned <- pseudo_native_eval(z, factor = 1)
  reference <- pseudo_native_reference(z, dose.factor = 1)
  expect_equal(aligned$grid.scores, reference$scores, tolerance = 1e-12)
})

test_that("ridge solves and intercept GCV match the explicit smoother", {
  z <- pseudo_native_fixture()
  for (intercept in c(FALSE, TRUE)) {
    actual <- pseudo_native_eval(z, intercept = intercept)
    expected <- pseudo_native_reference(z, intercept = intercept)
    expect_equal(drop(actual$GCV), expected$GCV, tolerance = 1e-10)
    expect_identical(actual$lambda.init, expected$lambda)
    expect_equal(drop(actual$W), expected$W, tolerance = 1e-10)
    expect_equal(actual$intercept, expected$alpha, tolerance = 1e-10)
    scalar <- z; scalar$lambda <- z$lambda[1]
    one <- pseudo_native_eval(scalar, intercept = intercept)
    expect_identical(one$lambda.init, scalar$lambda)
    if (intercept) expect_lt(abs(sum(actual$W)), 1e-12)
  }
  original <- pseudo_native_eval(z, intercept = TRUE)
  for (scale in c(.01, 3)) for (shift in c(0, 1e4)) {
    transformed <- z
    transformed$a <- scale * z$a + shift
    transformed$grid <- scale * z$grid + shift
    actual <- pseudo_native_eval(transformed, intercept = TRUE)
    expect_equal(drop(actual$grid.targets), scale * drop(original$grid.targets) + shift,
                 tolerance = 1e-12)
    expect_equal(actual$W, scale * original$W, tolerance = 1e-8)
    expect_equal(actual$intercept, scale * original$intercept + shift, tolerance = 1e-10)
    expect_equal(actual$GCV, scale^2 * original$GCV, tolerance = 1e-8)
    expect_identical(actual$lambda.init, original$lambda.init)
  }
})

test_that("constant-target intercepts and unusable ridge candidates have explicit behavior", {
  z <- pseudo_native_fixture(); z$grid <- .1; z$lambda <- c(.2, .03, .2)
  actual <- pseudo_native_eval(z, intercept = TRUE)
  expect_identical(actual$intercept, .1)
  expect_equal(drop(actual$W), rep(0, nrow(z$x)), tolerance = 0)
  expect_equal(drop(actual$GCV), rep(0, 3), tolerance = 0)
  expect_identical(actual$lambda.init, .2)
  # The large bandwidth makes the kernel exactly one in floating point.
  # A ridge below machine precision cannot make that represented matrix positive definite.
  z$bw <- 1e15; z$lambda <- c(1e-30, .1)
  expect_warning(actual <- pseudo_native_eval(z, intercept = TRUE), "Skipped unusable.*1e-30")
  expect_equal(actual$kernel, matrix(1, nrow(z$x), nrow(z$x)), tolerance = 0)
  expect_equal(drop(actual$ridge.valid), c(0, 1))
  expect_true(is.na(actual$GCV[1]))
  expect_identical(actual$lambda.init, .1)
  expect_identical(actual$intercept, .1)
  z$lambda <- 1e-30
  expect_error(pseudo_native_eval(z, intercept = TRUE), "No supplied.*ridge candidate")
})

test_that("optimizer convergence is gradient-only and accepted states are feasible", {
  for (p in c(3L, 6L)) {
    z <- pseudo_native_fixture(p)
    zero <- pseudo_native_audit(z, control = list(gtol = 1e-20))
    expect_identical(zero$stop.reason, "max_iterations")
    expect_equal(zero$itr, 0)
    expect_equal(zero$B, z$B, tolerance = 0)
    one <- pseudo_native_audit(z, 1L, control = list(gtol = 1e-20))
    expect_equal(one$itr, 1)
    expect_false(one$converge)
    expect_lt(norm(crossprod(one$B) - diag(2), "F"), 1e-10)
    expected <- pseudo_native_eval(z, one$B)
    expect_equal(one$objective.scaled, expected$objective.scaled, tolerance = 1e-12)
    expect_equal(one$gradient.norm, expected$gradient.norm, tolerance = 1e-12)
    stagnant <- pseudo_native_audit(z, 8L,
      control = list(gtol = 1e-20, btol = 10, ftol = 10))
    expect_equal(stagnant$itr, 5)
    expect_identical(stagnant$stop.reason, "stagnation")
    expect_false(stagnant$converge)
    stationary <- pseudo_native_audit(z, control = list(gtol = 1e6))
    expect_true(stationary$converge)
    expect_identical(stationary$stop.reason, "gradient_tolerance")
  }
})

test_that("thirty rejected trials roll back the objective, derivative and histories", {
  z <- pseudo_native_fixture()
  for (hook in c("trial", "gradient")) {
    arguments <- if (hook == "trial") list(fail.trial = 2L) else list(fail.gradient = 2L)
    previous <- pseudo_native_audit(z, 1L, control = list(gtol = 1e-20))
    failed <- do.call(pseudo_native_audit, c(list(z = z, maxitr = 3L,
      control = list(gtol = 1e-20)), arguments))
    expect_identical(failed$stop.reason, "line_search_failed")
    expect_false(failed$converge)
    expect_equal(failed$itr, 1)
    for (name in c("B", "W", "fn", "objective.scaled", "gradient.norm", "audit.history"))
      expect_equal(failed[[name]], previous[[name]], tolerance = 0)
    trials <- Filter(function(item) item$iteration == 2L, pseudo_native_events(failed, "trial"))
    expect_length(trials, 30L)
    expect_true(all(!vapply(trials, `[[`, logical(1), "accepted")))
    expect_length(pseudo_native_events(failed, "rollback"), 1L)
    for (trial in trials) expect_equal(trial$B.start, previous$B, tolerance = 0)
  }
  zero <- pseudo_native_audit(z, 2L, budget = 0L)
  expect_equal(zero$B, z$B, tolerance = 0)
  expect_equal(zero$itr, 0)
  expect_identical(zero$stop.reason, "line_search_failed")
})

test_that("nonmonotone reference updates and boundary gamma values are literal", {
  z <- pseudo_native_fixture()
  for (gamma in c(0, 1)) {
    fit <- pseudo_native_audit(z, 12L,
      control = list(gamma = gamma, gtol = 1e-20, tau = 30))
    initial <- pseudo_native_events(fit, "initialize")[[1]]
    values <- initial$objective
    for (commit in pseudo_native_events(fit, "commit")) {
      values <- c(values, commit$objective)
      expect_equal(commit$reference,
        if (gamma == 1) mean(values) else tail(values, 1), tolerance = 1e-12)
      expect_equal(commit$reference.weight, if (gamma == 1) length(values) else 1)
      expect_lt(norm(crossprod(commit$B) - diag(2), "F"), 1e-10)
    }
    for (trial in pseudo_native_events(fit, "trial")) if (trial$accepted) {
      expect_gt(trial$slope, 0)
      expect_lte(trial$objective, trial$reference - 1e-4 * trial$step * trial$slope)
    }
  }
})

test_that("invalid native geometry and controls fail clearly", {
  z <- pseudo_native_fixture()
  invalid <- z; invalid$B[, 2] <- invalid$B[, 1]
  expect_error(pseudo_native_audit(invalid), "orthonormal")
  invalid <- z; invalid$x[,] <- 1
  expect_error(pseudo_native_audit(invalid), "Projected covariates")
  invalid <- z; invalid$a[] <- .1
  expect_error(pseudo_native_audit(invalid), "Observed doses")
  for (control in list(list(gamma = -1), list(gamma = 1.1), list(eta = 1),
       list(rho = 0), list(tau = Inf), list(btol = 0), list(gtol = NA_real_)))
    expect_error(pseudo_native_audit(z, control = control), "Invalid pseudo-direct optimizer controls")
  for (lambda in list(0, -1, NA_real_, Inf, numeric())) {
    invalid <- z; invalid$lambda <- lambda
    expect_error(pseudo_native_eval(invalid), "lambda.*positive finite")
  }
})
test_that("large finite ridge values retain intercept GCV degrees of freedom", {
  set.seed(38191)
  n <- 17L
  x <- matrix(rnorm(n * 3), n, 3)
  a <- runif(n, 0, 4)
  r <- cos(a) + .7 * x[, 1] - .2 * x[, 2]^2
  evaluated <- orthoDr:::pdose_semi_evaluate(matrix(c(1, 0, 0), 3, 1),
    x, r, a, .7, seq(0, 4, length.out = 7), c(100, 1e100, 1e200), TRUE)
  target <- drop(evaluated$grid.targets)
  intercept_only <- n * sum((target - mean(target))^2) / (n - 1)^2
  expect_equal(drop(evaluated$GCV)[2:3], rep(intercept_only, 2), tolerance = 1e-12)
  expect_identical(evaluated$lambda.init, 100)
})

