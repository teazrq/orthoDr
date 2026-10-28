# Independent equations for the centered policy. These tests do not depend on
# helpers declared in other test files or on the native kernel implementation.
intercept_fixture <- function(d = 2L, repeated = FALSE) {
  n <- 12L
  set.seed(872L + d)
  x <- matrix(rnorm(3 * n), n, 3)
  x[, 2] <- .7 * x[, 1] + 1.4 * x[, 2]
  if (repeated) x[n, ] <- x[1, ]
  B <- qr.Q(qr(matrix(sin(seq_len(3 * d)), 3, d)))
  a <- seq(.1, 2.3, length.out = n)
  r <- -.4 + .3 * x[, 1] - (a - 1.2 - .6 * x[, 2])^2
  list(x = x, B = B, a = a, r = r, bw = .8,
       w = .13 + .09 * cos(seq_len(n)), alpha = 1.1)
}

intercept_reference <- function(z, B = z$B, w = z$w, alpha = z$alpha,
                                penalty = .17) {
  n <- nrow(z$x)
  u <- z$x %*% B
  s2 <- apply(u, 2, var)
  logk <- matrix(0, n, n)
  for (k in seq_len(ncol(B)))
    logk <- logk - outer(u[, k], u[, k], "-")^2 / (2 * z$bw^2 * s2[k])
  K <- exp(logk)
  C <- sweep(K, 2, colMeans(K), "-")
  dose <- drop(alpha + C %*% w)
  h2 <- var(z$a) * z$bw^2
  # Cancel the common dose^2 term before normalizing the reward weights.
  logits <- logk + outer(2 * z$a / h2, dose) - z$a^2 / h2
  P <- exp(sweep(logits, 2, apply(logits, 2, max), "-"))
  P <- sweep(P, 2, colSums(P), "/")
  mean.r <- drop(crossprod(P, z$r))
  g <- q <- numeric(n)
  for (i in seq_len(n)) {
    mean.a <- sum(P[, i] * z$a)
    g[i] <- -2 * sum(P[, i] * (z$r - mean.r[i]) * z$a) / (n * h2)
    q[i] <- -4 * sum(P[, i] * (z$r - mean.r[i]) * (z$a - mean.a)^2) / (n * h2^2)
  }
  F <- cbind(1, C)
  gradient <- drop(crossprod(F, g)) + c(0, 2 * penalty * w)
  hessian <- crossprod(F, F * q) + diag(c(0, rep(2 * penalty, n)))
  sensitivity <- -P * sweep(matrix(z$r, n, n), 2, mean.r, "-") / n +
    K * outer(w, g - mean(g))
  gB <- matrix(0, nrow(B), ncol(B))
  covariance <- cov(z$x)
  for (k in seq_len(ncol(B))) for (i in seq_len(n)) for (j in seq_len(n)) {
    delta <- u[j, k] - u[i, k]
    pair <- -delta * (z$x[j, ] - z$x[i, ]) / (z$bw^2 * s2[k]) +
      delta^2 * drop(covariance %*% B[, k]) / (z$bw^2 * s2[k]^2)
    gB[, k] <- gB[, k] + sensitivity[j, i] * pair
  }
  list(fn = -mean(mean.r) + penalty * sum(w^2), value = mean(mean.r),
       penalty = penalty * sum(w^2), kernel = K, features = C,
       dose = dose, gradient = gradient, hessian = hessian, gB = gB)
}

intercept_native <- function(z, B = z$B, w = z$w, alpha = z$alpha,
                             penalty = .17, bounds = range(z$a)) {
  orthoDr:::pdose_direct_evaluate(B, z$x, z$a, z$r, z$bw, w, penalty,
                                  0, TRUE, bounds, alpha)
}

intercept_audit <- function(z, maxitr = 0L, lambda = .1, penalty = .17,
                            bounds = range(z$a), w = NULL, alpha = NULL,
                            grid = seq(min(z$a), max(z$a), length.out = 7L),
                            B.budget = 30L, control = list()) {
  ctl <- modifyList(list(rho = 1e-4, eta = .2, tau = 1,
                         btol = 1e-10, ftol = 1e-10, gtol = 1e-6), control)
  orthoDr:::pdose_direct_audit_solver(z$B, z$x, z$a,
    exp(-outer(z$a, grid, "-")^2 / (sd(z$a) * z$bw)^2), grid,
    z$r, lambda, z$bw, ctl$rho, ctl$eta, .85, ctl$tau, 1e-6,
    ctl$btol, ctl$ftol, ctl$gtol, maxitr, 0L, 1L, penalty,
    TRUE, B.budget, w, TRUE, bounds, alpha)
}

intercept_events <- function(fit, event, block = NULL) {
  Filter(function(item) identical(item$event, event) &&
           (is.null(block) || identical(item$block, block)), fit$trace)
}

intercept_expect_state <- function(fit, z, penalty = .17, gtol = 1e-6) {
  ref <- intercept_reference(z, fit$B, drop(fit$W), fit$intercept, penalty)
  expect_equal(fit$fn, ref$fn, tolerance = 1e-10)
  expect_equal(fit$value, ref$value, tolerance = 1e-10)
  expect_equal(fit$penalty, ref$penalty, tolerance = 1e-10)
  expect_equal(drop(fit$kernel.mean), colMeans(ref$kernel), tolerance = 1e-12)
  ga <- ref$gradient[1]
  bounds <- fit$intercept.bounds
  residual <- fit$intercept - min(bounds[2], max(bounds[1], fit$intercept - ga))
  gw <- sqrt(sum(ref$gradient[-1]^2))
  tangent <- ref$gB - fit$B %*% ((crossprod(fit$B, ref$gB) + crossprod(ref$gB, fit$B)) / 2)
  expect_equal(fit$gradient.intercept, ga, tolerance = 1e-10)
  expect_equal(unname(fit$gradient.norm),
               c(norm(tangent, "F"), gw, abs(residual), sqrt(gw^2 + residual^2)),
               tolerance = 1e-9)
  expect_identical(fit$converge,
    fit$gradient.norm[["B"]] <= gtol && fit$gradient.norm[["policy"]] <= gtol)
  expect_equal(mean(ref$dose), fit$intercept, tolerance = 1e-12)
  expect_gte(fit$intercept, bounds[1])
  expect_lte(fit$intercept, bounds[2])
}

test_that("centered policy gradients, joint Hessian and moving means match independent equations", {
  for (d in 1:2) {
    z <- intercept_fixture(d)
    B <- z$B + matrix(seq_along(z$B), nrow(z$B)) / 140
    native <- intercept_native(z, B = B)
    ref <- intercept_reference(z, B = B)
    expect_true(native$valid)
    for (name in c("fn", "value", "penalty", "kernel", "features", "dose"))
      expect_equal(native[[name]], ref[[name]], tolerance = 1e-11, ignore_attr = TRUE)
    expect_equal(c(native$g.intercept, drop(native$gW)), ref$gradient, tolerance = 1e-10)
    expect_equal(native$hessian.policy, ref$hessian, tolerance = 1e-10)
    expect_equal(native$hessian.W, ref$hessian[-1, -1], tolerance = 1e-10)
    expect_equal(native$gB, ref$gB, tolerance = 1e-10)
    expect_equal(mean(native$dose), z$alpha, tolerance = 1e-12)
    expect_gt(abs(sum(z$w)), .1)
    theta <- c(z$alpha, z$w)
    eps <- 2e-6
    gfd <- numeric(length(theta))
    Hfd <- matrix(0, length(theta), length(theta))
    for (j in seq_along(theta)) {
      delta <- numeric(length(theta)); delta[j] <- eps
      plus <- theta + delta; minus <- theta - delta
      zp <- intercept_reference(z, B, plus[-1], plus[1])
      zm <- intercept_reference(z, B, minus[-1], minus[1])
      gfd[j] <- (zp$fn - zm$fn) / (2 * eps)
      Hfd[, j] <- (zp$gradient - zm$gradient) / (2 * eps)
    }
    expect_equal(ref$gradient, gfd, tolerance = 1e-7)
    expect_equal(ref$hessian, Hfd, tolerance = 1e-7)
    D <- matrix(cos(seq_along(B)), nrow(B))
    derivative <- (intercept_reference(z, B + eps * D)$fn -
                    intercept_reference(z, B - eps * D)$fn) / (2 * eps)
    expect_equal(sum(native$gB * D), derivative, tolerance = 1e-7)
  }
})

test_that("feasible Newton proposals recheck clipped slopes and use continuous residuals", {
  H <- matrix(c(5, 2, 2, 1), 2, 2)
  g <- c(10, 1)
  raw <- -drop(solve(H, g))
  clipped <- c(max(-.5, min(.5, raw[1])), raw[2])
  expect_equal(sum(g * raw), -65)
  expect_equal(sum(g * clipped), 10)
  native <- orthoDr:::pdose_direct_policy_direction(g, H, 0, c(-.5, .5))
  expect_identical(native$direction.type, "gradient")
  expect_equal(drop(native$direction), c(-.5, -1))
  expect_equal(native$slope, -6)
  for (alpha in c(-.5, .5)) {
    sign <- if (alpha < 0) 1 else -1
    active <- orthoDr:::pdose_direct_policy_direction(c(sign, 2), diag(c(0, 4)), alpha, c(-.5, .5))
    expect_identical(active$direction.type, "newton.fixed.intercept")
    expect_equal(drop(active$direction), c(0, -.5))
    expect_equal(active$residual.intercept, 0)
  }
  near <- orthoDr:::pdose_direct_policy_direction(c(1, 0), diag(2),
                                                  -1 + 1.5e-9, c(-1, 1))
  expect_lt(abs(near$residual.intercept), 2e-9)
  expect_equal(drop(near$direction), c(-1.5e-9, 0), tolerance = 1e-7)
})

test_that("intercept KRR and GCV include stable contrast degrees of freedom", {
  for (repeated in c(FALSE, TRUE)) {
    z <- intercept_fixture(repeated = repeated)
    candidates <- c(.03, .2, .03, 1)
    fit <- intercept_audit(z, lambda = candidates)
    initial <- intercept_events(fit, "initialize")[[1]]
    t <- drop(initial$grid.targets)
    expect_gt(length(unique(t)), 1L)
    n <- length(t)
    K <- intercept_reference(z)$kernel
    H <- diag(n) - matrix(1 / n, n, n)
    G <- H %*% K %*% H
    Z <- qr.Q(qr(cbind(rep(1, n), diag(n)[, -1])))[, -1]
    scores <- vapply(candidates, function(lambda) {
      inverse <- solve(crossprod(Z, K %*% Z) + lambda * diag(n - 1))
      n * sum((inverse %*% crossprod(Z, t))^2) / sum(diag(inverse))^2
    }, numeric(1))
    expect_equal(drop(initial$GCV), scores, tolerance = 1e-10)
    selected <- candidates[which.min(scores)]
    expect_equal(fit$lambda.init, selected)
    expect_equal(drop(fit$W), drop(solve(G + selected * diag(n), H %*% t)), tolerance = 1e-10)
    expect_equal(fit$intercept, mean(t), tolerance = 1e-12)
    expect_equal(sum(fit$W), 0, tolerance = 1e-12)
    for (lambda in candidates) {
      S <- matrix(1 / n, n, n) + G %*% solve(G + lambda * diag(n))
      score <- mean((t - S %*% t)^2) / (1 - sum(diag(S)) / n)^2
      expect_equal(score, scores[match(lambda, candidates)], tolerance = 1e-9)
    }
    # Tiny lambda on a diagonal kernel checks cancellation-free residual df.
    diagonal <- list(x = matrix(c(-1, 0, 1)), B = matrix(1), a = c(0, 1, 2),
                     r = c(1, 2, 3), bw = .001)
    stable <- intercept_audit(diagonal, lambda = c(1e-20, .1), grid = c(0, 1, 2))
    expect_true(all(is.finite(intercept_events(stable, "initialize")[[1]]$GCV)))
  }
})

test_that("constant targets are exact and custom active bounds require scalar initialization", {
  z <- intercept_fixture()
  for (bounds in list(c(.1, 2), c(-1, .1), c(0, 1))) {
    fit <- intercept_audit(z, lambda = c(.2, .03), grid = .1, bounds = bounds)
    expect_identical(fit$lambda.init, .2)
    expect_identical(fit$intercept, .1)
    expect_equal(drop(fit$W), rep(0, nrow(z$x)), tolerance = 0)
    expect_equal(drop(intercept_events(fit, "initialize")[[1]]$GCV), c(0, 0))
  }
  expect_error(intercept_audit(z, lambda = c(.1, .2), bounds = c(-2, -1)),
               "Active intercept.bounds.*scalar lambda")
  expect_error(intercept_audit(z, lambda = c(.1, .2), bounds = c(-2, -1), grid = .1),
               "Active intercept.bounds.*scalar lambda")
  constrained <- intercept_audit(z, bounds = c(-2, -1))
  initial <- intercept_events(constrained, "initialize")[[1]]
  n <- nrow(z$x); H <- diag(n) - matrix(1 / n, n, n)
  G <- H %*% intercept_reference(z)$kernel %*% H
  expect_identical(constrained$intercept, -1)
  expect_equal(drop(constrained$W), drop(solve(G + .1 * diag(n), H %*% initial$grid.targets)),
               tolerance = 1e-10)
  intercept_expect_state(constrained, z)
})

test_that("bounded stationarity permits an outward raw gradient for monotone rewards", {
  z <- list(x = matrix(c(-1, -1, 1, 1)), B = matrix(1), a = c(0, 2, 0, 2),
            r = c(0, 2, 0, 2), bw = 1)
  fit <- intercept_audit(z, w = rep(0, 4), alpha = 2)
  expect_identical(fit$stop.reason, "gradient_tolerance")
  expect_true(fit$converge)
  expect_lt(fit$gradient.intercept, -.1)
  expect_equal(fit$gradient.norm[["intercept"]], 0)
  intercept_expect_state(fit, z)
  near <- intercept_audit(z, w = rep(0, 4), alpha = 2 - 1.5e-9,
                          control = list(gtol = 2e-9))
  expect_true(near$converge)
  expect_lt(near$gradient.intercept, -.1)
  expect_lt(near$gradient.norm[["policy"]], 2e-9)
})

test_that("bounded policy and direction trials are monotone and roll back atomically", {
  z <- intercept_fixture()
  for (tau in c(1e-4, 1000)) {
    fit <- intercept_audit(z, maxitr = 2, w = z$w, alpha = z$alpha,
                            control = list(tau = tau, gtol = 1e-12))
    trials <- intercept_events(fit, "trial")
    expect_gt(length(trials), 0L)
    accepted <- Filter(function(item) isTRUE(item$accepted), trials)
    expect_gt(length(accepted), 0L)
    for (item in trials) {
      start <- intercept_reference(z, item$B.start, drop(item$W.start), item$intercept.start)
      expect_equal(item$fn.start, start$fn, tolerance = 1e-10)
      expect_gte(item$intercept, min(z$a))
      expect_lte(item$intercept, max(z$a))
      if (item$finite) {
        ref <- intercept_reference(z, item$B, drop(item$W), item$intercept)
        expect_equal(item$fn, ref$fn, tolerance = 1e-10)
        expect_equal(drop(item$kernel.mean), colMeans(ref$kernel), tolerance = 1e-12)
      }
      if (item$block == "policy") {
        expect_lte(item$step, 1)
        actual <- sum(start$gradient * c(item$direction.intercept, item$direction.W))
        expect_equal(item$slope, actual, tolerance = 1e-10)
      } else expect_equal(item$intercept, item$intercept.start)
      if (item$accepted) {
        expect_lt(item$slope, 0)
        expect_lte(item$fn, item$fn.start + 1e-4 * item$step * item$slope + 1e-12)
      }
    }
    expect_length(intercept_events(fit, "commit"), fit$itr)
    intercept_expect_state(fit, z, gtol = 1e-12)
  }
  rollback <- intercept_audit(z, maxitr = 1, w = z$w, alpha = z$alpha,
    B.budget = 0L, control = list(tau = 1e-4, gtol = 1e-12))
  accepted <- Filter(function(item) item$accepted, intercept_events(rollback, "trial", "policy"))
  expect_gt(length(accepted), 0L)
  expect_true(any(vapply(accepted, function(item) item$intercept != z$alpha, logical(1))))
  expect_identical(rollback$stop.reason, "B_line_search_failed")
  expect_length(intercept_events(rollback, "rollback"), 1L)
  expect_equal(rollback$B, z$B, tolerance = 0)
  expect_equal(drop(rollback$W), z$w, tolerance = 0)
  expect_identical(rollback$intercept, z$alpha)
  expect_equal(rollback$itr, 0)
  intercept_expect_state(rollback, z, gtol = 1e-12)
})

test_that("constant rewards, finite budgets and stagnation retain fresh policy diagnostics", {
  z <- intercept_fixture()
  flat <- z; flat$r[] <- -3
  fit <- intercept_audit(flat, maxitr = 1, penalty = .5, w = z$w, alpha = z$alpha,
                          control = list(tau = 1, gtol = 1e-10))
  expect_true(fit$converge)
  expect_identical(fit$intercept, z$alpha)
  expect_lt(sqrt(sum(fit$W^2)), 1e-10)
  intercept_expect_state(fit, flat, penalty = .5, gtol = 1e-10)
  initial <- intercept_audit(z, w = z$w, alpha = z$alpha)
  expect_identical(initial$stop.reason, "iteration_limit")
  intercept_expect_state(initial, z)
  small <- intercept_audit(z, maxitr = 8, w = z$w, alpha = z$alpha,
    # Loose change thresholds deliberately exercise stagnation. A tiny
    # initial tau no longer implies tiny later steps after spectral scaling.
    control = list(tau = 1e-12, gtol = 1e-20, btol = 10, ftol = 10))
  expect_equal(small$itr, 5)
  expect_identical(small$stop.reason, "stagnation")
  expect_false(small$converge)
  intercept_expect_state(small, z, gtol = 1e-20)
})

test_that("dose and reward affine transformations preserve the stated fixed-state identities", {
  z <- intercept_fixture()
  original <- intercept_native(z)
  transformed <- z
  transformed$a <- 3 * z$a + 8
  transformed$w <- 3 * z$w
  transformed$alpha <- 3 * z$alpha + 8
  fit <- intercept_native(transformed, penalty = .17 / 9)
  expect_equal(fit$fn, original$fn, tolerance = 1e-11)
  expect_equal(fit$dose, 3 * original$dose + 8, tolerance = 1e-11)
  initial <- intercept_audit(z)
  shifted <- intercept_audit(transformed)
  expect_equal(shifted$intercept, 3 * initial$intercept + 8, tolerance = 1e-10)
  expect_equal(shifted$W, 3 * initial$W, tolerance = 1e-10)
  transformed <- z
  transformed$r <- 2 * z$r - 4
  reward <- intercept_native(transformed, penalty = 2 * .17)
  expect_equal(reward$fn, 2 * original$fn + 4, tolerance = 1e-11)
})

test_that("native intercept inputs reject zero penalty and invalid custom initial states", {
  z <- intercept_fixture()
  expect_error(intercept_audit(z, penalty = 0), "lambda.reward.*positive")
  expect_false(intercept_native(z, penalty = 0)$valid)
  for (bounds in list(c(1, 1), c(2, 1), c(NA, 1), c(-Inf, 1), 1))
    expect_error(intercept_audit(z, bounds = bounds), "intercept.bounds")
  for (alpha in list(NA_real_, Inf, -10, c(1, 2)))
    expect_error(intercept_audit(z, w = z$w, alpha = alpha), "alpha_initial")
  legacy <- orthoDr:::pdose_direct_evaluate(z$B, z$x, z$a, z$r, z$bw, z$w, 0)
  expect_true(legacy$valid)
  expect_named(legacy$gradient.norm, c("B", "W"))
})
