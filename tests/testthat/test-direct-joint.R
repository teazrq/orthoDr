# Equation-based references deliberately use no native kernel/objective helpers.
direct_fixture <- function(p = 4L, d = 2L, n = 14L) {
  set.seed(318L + p + d)
  x <- matrix(rnorm(n * p), n, p)
  x[, 2] <- .7 * x[, 1] + 1.8 * x[, 2]
  x <- sweep(x, 2, seq(.6, 2, length.out = p), "*")
  B <- qr.Q(qr(matrix(sin(seq_len(p * d)) + cos(2 * seq_len(p * d)), p, d)))
  a <- seq(.2, 2.8, length.out = n) + .05 * x[, 1]
  r <- -.7 + .4 * x[, 1] - (a - 1 - .2 * x[, 2])^2
  list(x = x, a = a, r = r, B = B, w = .1 + .17 * sin(seq_len(n)), bw = .65)
}

direct_reference <- function(B, w, x, a, r, bw, lambda.reward = .1) {
  n <- nrow(x)
  u <- x %*% B
  scales <- vapply(seq_len(ncol(B)), function(l) sd(u[, l]), numeric(1))
  if (any(!is.finite(scales) | scales <= 0)) stop("invalid projected scale")
  h2 <- (bw * sd(a))^2
  logk <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) for (l in seq_len(ncol(B)))
    logk[i, j] <- logk[i, j] - (u[i, l] - u[j, l])^2 / (2 * bw^2 * scales[l]^2)
  kernel <- exp(logk)
  dose <- drop(kernel %*% w)
  prob <- matrix(0, n, n)
  mean.reward <- numeric(n)
  for (i in seq_len(n)) {
    lw <- numeric(n)
    for (j in seq_len(n)) lw[j] <- logk[j, i] - (a[j] - dose[i])^2 / h2
    relative <- exp(lw - max(lw))
    prob[, i] <- relative / sum(relative)
    mean.reward[i] <- sum(prob[, i] * r)
  }
  value <- mean(mean.reward)
  penalty <- lambda.reward * sum(w^2)
  list(fn = -value + penalty, value = value, penalty = penalty,
       kernel = kernel, logkernel = logk, dose = dose, scales = scales, prob = prob,
       mean.reward = mean.reward, u = u, h2 = h2)
}

direct_reference_w <- function(B, w, x, a, r, bw, lambda.reward = .1) {
  z <- direct_reference(B, w, x, a, r, bw, lambda.reward)
  gdose <- numeric(nrow(x))
  for (i in seq_len(nrow(x))) for (j in seq_len(nrow(x)))
    gdose[i] <- gdose[i] - 2 / (nrow(x) * z$h2) *
      z$prob[j, i] * (r[j] - z$mean.reward[i]) * (a[j] - z$dose[i])
  drop(crossprod(z$kernel, gdose)) + 2 * lambda.reward * w
}

direct_reference_hessian <- function(B, w, x, a, r, bw, lambda.reward = .1) {
  z <- direct_reference(B, w, x, a, r, bw, lambda.reward)
  curvature <- numeric(nrow(x))
  for (i in seq_len(nrow(x))) {
    mean.a <- sum(z$prob[, i] * a)
    covariance <- 0
    for (j in seq_len(nrow(x)))
      covariance <- covariance + z$prob[j, i] * (r[j] - z$mean.reward[i]) *
        ((a[j] - mean.a) / sqrt(z$h2))^2
    curvature[i] <- -4 * covariance / (nrow(x) * z$h2)
  }
  crossprod(z$kernel, curvature * z$kernel) + diag(2 * lambda.reward, nrow(x))
}

# Differentiate the equations along an arbitrary direction, including the
# moving sample centres and sample standard deviations.
direct_reference_direction <- function(B, w, D, x, a, r, bw) {
  z <- direct_reference(B, w, x, a, r, bw, 0)
  du <- x %*% D
  ds <- vapply(seq_len(ncol(B)), function(l) cov(z$u[, l], du[, l]) / z$scales[l], numeric(1))
  dlk <- matrix(0, nrow(x), nrow(x))
  for (i in seq_len(nrow(x))) for (j in seq_len(nrow(x))) for (l in seq_len(ncol(B))) {
    delta <- z$u[i, l] - z$u[j, l]
    ddelta <- du[i, l] - du[j, l]
    dlk[i, j] <- dlk[i, j] - delta * ddelta / (bw^2 * z$scales[l]^2) +
      delta^2 * ds[l] / (bw^2 * z$scales[l]^3)
  }
  ddose <- drop((z$kernel * dlk) %*% w)
  answer <- 0
  for (i in seq_len(nrow(x))) for (j in seq_len(nrow(x)))
    answer <- answer - z$prob[j, i] * (r[j] - z$mean.reward[i]) / nrow(x) *
      (dlk[j, i] + 2 * (a[j] - z$dose[i]) * ddose[i] / z$h2)
  answer
}

direct_reference_B <- function(B, w, x, a, r, bw) {
  out <- matrix(0, nrow(B), ncol(B))
  for (k in seq_along(B)) {
    direction <- matrix(0, nrow(B), ncol(B))
    direction[k] <- 1
    out[k] <- direct_reference_direction(B, w, direction, x, a, r, bw)
  }
  out
}

direct_reference_predict <- function(B, w, x, newx, bw) {
  train <- x %*% B
  test <- newx %*% B
  scales <- vapply(seq_len(ncol(B)), function(l) sd(train[, l]), numeric(1))
  out <- numeric(nrow(newx))
  for (i in seq_len(nrow(newx))) for (j in seq_len(nrow(x))) {
    exponent <- 0
    for (l in seq_len(ncol(B)))
      exponent <- exponent - (test[i, l] - train[j, l])^2 / (2 * bw^2 * scales[l]^2)
    out[i] <- out[i] + exp(exponent) * w[j]
  }
  out
}

direct_native <- function(z, B = z$B, w = z$w, r = z$r, a = z$a,
                          penalty = .1, step = 0) {
  orthoDr:::pdose_direct_evaluate(B, z$x, a, r, z$bw, w, penalty, step)
}

direct_expect_close <- function(actual, reference) {
  expect_lte(max(abs(as.vector(actual) - as.vector(reference))),
             1e-10 * max(1, abs(reference)))
}

direct_expect_differences <- function(fun, reference) {
  steps <- c(1e-4, 1e-5, 1e-6)
  differences <- vapply(steps, function(h) (fun(h) - fun(-h)) / (2 * h), numeric(1))
  agrees <- abs(differences - reference) <= 1e-7 + 1e-5 * abs(reference)
  expect_true(any(head(agrees, -1L) & tail(agrees, -1L)),
              info = paste("central differences:", paste(signif(differences, 9), collapse = ", ")))
}

test_that("direct objective and both gradients agree with independent equations", {
  for (p in c(2L, 4L, 6L)) for (penalty in c(0, .17)) {
    z <- direct_fixture(p)
    # Arbitrary coefficients and a nonorthogonal evaluation state test more
    # than a ridge start or a constrained stationary point.
    B <- z$B + matrix(seq_len(length(z$B)), nrow(z$B)) / 130
    native <- direct_native(z, B = B, penalty = penalty)
    reference <- direct_reference(B, z$w, z$x, z$a, z$r, z$bw, penalty)
    expect_true(native$valid)
    for (component in c("fn", "value", "penalty", "kernel", "dose", "scales"))
      direct_expect_close(native[[component]], reference[[component]])
    gw <- direct_reference_w(B, z$w, z$x, z$a, z$r, z$bw, penalty)
    direct_expect_close(native$gW, gw)
    for (k in 1:2) {
      D <- matrix(sin(k * seq_along(B)), nrow(B), ncol(B))
      v <- cos(k * seq_along(z$w))
      derivative.B <- direct_reference_direction(B, z$w, D, z$x, z$a, z$r, z$bw)
      derivative.W <- sum(gw * v)
      expect_lte(abs(sum(native$gB * D) - derivative.B), 1e-7 + 1e-5 * abs(derivative.B))
      direct_expect_differences(function(h)
        direct_reference(B + h * D, z$w, z$x, z$a, z$r, z$bw, penalty)$fn, derivative.B)
      direct_expect_differences(function(h)
        direct_reference(B, z$w + h * v, z$x, z$a, z$r, z$bw, penalty)$fn, derivative.W)
      direct_expect_differences(function(h)
        direct_reference(B + h * D, z$w + h * v, z$x, z$a, z$r, z$bw, penalty)$fn,
        derivative.B + derivative.W)
    }
  }
})

test_that("full Stiefel residual and both Cayley branches use their actual slope", {
  for (p in c(2L, 4L, 6L)) {
    z <- direct_fixture(p)
    if (p == 2L) z$B <- diag(2)
    native <- direct_native(z, step = .003)
    G <- direct_reference_B(z$B, z$w, z$x, z$a, z$r, z$bw)
    BtG <- crossprod(z$B, G)
    projected <- G - z$B %*% ((BtG + t(BtG)) / 2)
    velocity <- -(G - z$B %*% crossprod(G, z$B))
    direct_expect_close(native$gradient.norm[["B"]], norm(projected, "F"))
    direct_expect_close(native$cayley.velocity, velocity)
    direct_expect_close(native$cayley.slope, sum(G * velocity))
    expect_equal(crossprod(z$B, projected) + crossprod(projected, z$B), matrix(0, 2, 2), tolerance = 1e-10)
    skew <- G %*% t(z$B) - z$B %*% t(G)
    curve <- function(t) solve(diag(p) + t / 2 * skew, (diag(p) - t / 2 * skew) %*% z$B)
    direct_expect_close(native$cayley.B, curve(.003))
    expect_equal(crossprod(native$cayley.B), diag(2), tolerance = 1e-10)
    direct_expect_differences(function(h)
      direct_reference(curve(h), z$w, z$x, z$a, z$r, z$bw)$fn, native$cayley.slope)
    if (p == 2L) {
      horizontal <- G - z$B %*% crossprod(z$B, G)
      expect_lt(norm(horizontal, "F"), 1e-12)
      expect_gt(norm(projected, "F"), 1e-4)
      expect_gt(abs(native$cayley.slope + sum(projected^2)), 1e-6)
    }
  }
})

test_that("direct coefficient Hessians agree with full independent gradient differences", {
  for (p in c(2L, 4L, 6L)) for (penalty in c(0, .17)) {
    z <- direct_fixture(p)
    B <- z$B + matrix(seq_along(z$B), nrow(z$B)) / 130
    native <- direct_native(z, B = B, penalty = penalty)
    H <- direct_reference_hessian(B, z$w, z$x, z$a, z$r, z$bw, penalty)
    expect_true(native$hessian.valid)
    direct_expect_close(native$hessian.W, H)
    expect_equal(H, t(H), tolerance = 1e-12)
    agreement <- logical(3)
    steps <- c(1e-4, 1e-5, 1e-6)
    for (k in seq_along(steps)) {
      differences <- matrix(NA_real_, length(z$w), length(z$w))
      for (j in seq_along(z$w)) {
        perturbation <- z$w * 0
        perturbation[j] <- steps[k]
        differences[, j] <-
          (direct_reference_w(B, z$w + perturbation, z$x, z$a, z$r, z$bw, penalty) -
             direct_reference_w(B, z$w - perturbation, z$x, z$a, z$r, z$bw, penalty)) / (2 * steps[k])
      }
      agreement[k] <- all(abs(H - differences) <= 1e-7 + 1e-5 * abs(H))
    }
    expect_true(any(head(agreement, -1L) & tail(agreement, -1L)))
    flat <- direct_native(z, B = B, r = rep(-3, nrow(z$x)), penalty = penalty)
    direct_expect_close(flat$hessian.W, diag(2 * penalty, nrow(z$x)))
    scaled <- direct_native(z, B = B, r = 3 * z$r + 7, penalty = 3 * penalty)
    direct_expect_close(scaled$hessian.W, 3 * H)
    units <- direct_native(z, B = B, a = 5 * z$a, w = 5 * z$w, penalty = penalty / 25)
    direct_expect_close(25 * units$hessian.W, H)
  }
})

test_that("direct log weights survive raw kernel and dose-weight underflow", {
  z <- list(x = matrix(c(-1, 0, 1)), a = c(0, 1, 2), r = c(-4, 2, 7),
            B = matrix(1), w = rep(1, 3), bw = .01)
  native <- direct_native(z, penalty = 0)
  reference <- direct_reference(z$B, z$w, z$x, z$a, z$r, z$bw, 0)
  raw <- reference$kernel * exp(-(outer(z$a, reference$dose, "-") / (sd(z$a) * z$bw))^2)
  expect_true(any(colSums(raw) == 0))
  expect_true(native$valid)
  direct_expect_close(native$fn, reference$fn)
  direct_expect_close(native$value, 2)
  expect_true(all(is.finite(native$gB)))
  expect_true(all(is.finite(native$gW)))
  expect_false(direct_native(z, B = matrix(0))$valid)
  expect_false(direct_native(z, w = rep(1e308, 3), penalty = 1)$valid)
})

test_that("direct penalty and reward/dose scale identities hold at a fixed state", {
  z <- direct_fixture()
  lambda <- .23
  base <- direct_native(z, penalty = lambda)
  transformed <- direct_native(z, r = 3 * z$r + 7, penalty = 3 * lambda)
  direct_expect_close(transformed$fn, 3 * base$fn - 7)
  direct_expect_close(transformed$gB, 3 * base$gB)
  direct_expect_close(transformed$gW, 3 * base$gW)
  units <- direct_native(z, a = 5 * z$a, w = 5 * z$w, penalty = lambda / 25)
  direct_expect_close(units$fn, base$fn)
  direct_expect_close(units$gB, base$gB)
  direct_expect_close(5 * units$gW, base$gW)
  flat <- direct_native(z, r = rep(-3, nrow(z$x)), penalty = lambda)
  direct_expect_close(flat$value, -3)
  direct_expect_close(flat$fn, 3 + lambda * sum(z$w^2))
  direct_expect_close(flat$gW, 2 * lambda * z$w)
  direct_expect_close(flat$gB, matrix(0, nrow(z$B), ncol(z$B)))
  # Sign and column permutation preserve the dictionary and policy exactly.
  equivalent <- direct_native(z, B = z$B[, 2:1, drop = FALSE] %*% diag(c(-1, 1)), penalty = lambda)
  direct_expect_close(equivalent$fn, base$fn)
  direct_expect_close(equivalent$dose, base$dose)
  direct_expect_close(equivalent$gB, base$gB[, 2:1, drop = FALSE] %*% diag(c(-1, 1)))
  direct_expect_close(equivalent$gW, base$gW)
})

direct_initial_reference <- function(z, lambda = .1, K = 5L, r = z$r) {
  grid <- seq(min(z$a), max(z$a), length.out = K)
  state <- direct_reference(z$B, rep(0, nrow(z$x)), z$x, z$a, r, z$bw, 0)
  targets <- numeric(nrow(z$x))
  for (i in seq_len(nrow(z$x))) {
    means <- numeric(K)
    for (k in seq_len(K)) {
      lw <- state$logkernel[, i] - (z$a - grid[k])^2 / state$h2
      weights <- exp(lw - max(lw))
      means[k] <- sum(weights * r) / sum(weights)
    }
    targets[i] <- grid[which.max(means)]
  }
  scores <- vapply(lambda, function(candidate) {
    residual.operator <- diag(nrow(z$x)) - state$kernel %*%
      solve(state$kernel + diag(candidate, nrow(z$x)))
    nrow(z$x) * sum((residual.operator %*% targets)^2) / sum(diag(residual.operator))^2
  }, numeric(1))
  chosen <- lambda[which.min(scores)]
  list(W = drop(solve(state$kernel + diag(chosen, nrow(z$x)), targets)),
       lambda = chosen, GCV = scores, targets = targets, grid = grid)
}

direct_audit <- function(z, maxitr = 1L, lambda = .1, penalty = .1,
                         K = 5L, control = list(), initial.w = NULL,
                         B.budget = 30L, r = z$r, ncore = 1L) {
  defaults <- list(rho = 1e-4, eta = .2, gamma = .85, tau = .001,
                   epsilon = 1e-6, btol = 1e-6, ftol = 1e-6, gtol = 1e-6)
  defaults[names(control)] <- control
  grid <- seq(min(z$a), max(z$a), length.out = K)
  args <- c(list(B = z$B, X = z$x, A = z$a,
                 a_dist = exp(-(outer(z$a, grid, "-") / (sd(z$a) * z$bw))^2),
                 a_seq = grid, R = r, lambda = lambda, bw = z$bw), defaults,
            list(maxitr = maxitr, verbose = 0L, ncore = ncore,
                 lambda_reward = penalty, trace = TRUE,
                 B_trial_budget = B.budget, W_initial = initial.w))
  do.call(orthoDr:::pdose_direct_audit_solver, args)
}

direct_public <- function(z, ...) {
  orthoDr_pdose(z$x, z$a, z$r, ndr = ncol(z$B), B.initial = z$B,
                bw = z$bw, method = "direct", ncore = 1L, ...)
}

direct_events <- function(fit, event, block = NULL) {
  Filter(function(item) identical(item$event, event) &&
           (is.null(block) || identical(item$block, block)), fit$trace)
}

direct_expect_returned_state <- function(fit, z, penalty = .1, r = z$r, gtol = 1e-6) {
  ref <- direct_reference(fit$B, drop(fit$W), z$x, z$a, r, z$bw, penalty)
  for (component in c("fn", "value", "penalty")) direct_expect_close(fit[[component]], ref[[component]])
  G <- direct_reference_B(fit$B, drop(fit$W), z$x, z$a, r, z$bw)
  BtG <- crossprod(fit$B, G)
  residual.B <- norm(G - fit$B %*% ((BtG + t(BtG)) / 2), "F")
  residual.W <- sqrt(sum(direct_reference_w(fit$B, drop(fit$W), z$x, z$a, r, z$bw, penalty)^2))
  direct_expect_close(fit$gradient.norm[["B"]], residual.B)
  direct_expect_close(fit$gradient.norm[["W"]], residual.W)
  expect_identical(fit$converge, identical(fit$stop.reason, "gradient_tolerance"))
  expect_identical(fit$converge, residual.B <= gtol && residual.W <= gtol)
  expect_equal(crossprod(fit$B), diag(ncol(fit$B)), tolerance = 1e-10)
  expect_identical(fit$objective.version, "direct_joint_v1")
}

test_that("direct grid targets, GCV and raw ridge inflation occur once at initialization", {
  z <- direct_fixture()
  for (K in c(1L, 5L)) for (lambda in list(.1, c(.03, .1, .5), c(.5, .03, .1))) {
    ref <- direct_initial_reference(z, lambda, K)
    fit <- direct_audit(z, maxitr = 0, lambda = lambda, K = K)
    expect_equal(drop(fit$W), ref$W, tolerance = 1e-10)
    expect_equal(fit$lambda.init, ref$lambda)
    expect_equal(fit$lambda.reward, .1)
    initialization <- direct_events(fit, "initialize")
    expect_length(initialization, 1L)
    expect_equal(initialization[[1]]$initialization.count, 1L)
    expect_equal(drop(initialization[[1]]$grid.targets), ref$targets)
    if (length(lambda) == 1L) {
      expect_length(initialization[[1]]$GCV, 0L)
    } else {
      expect_equal(drop(initialization[[1]]$GCV), ref$GCV, tolerance = 1e-10)
    }
    expect_equal(fit$itr, 0L)
    expect_identical(fit$stop.reason, "iteration_limit")
    direct_expect_returned_state(fit, z)
  }
  # Every GCV score is exactly zero when K=1 and its target is zero.
  tie <- z
  tie$a <- tie$a - min(tie$a)
  tied <- direct_audit(tie, lambda = c(.7, .03, .1), K = 1L, maxitr = 0)
  expect_equal(tied$lambda.init, .7)
  expect_equal(drop(tied$W), rep(0, nrow(z$x)))
  expect_equal(drop(direct_events(tied, "initialize")[[1]]$GCV), c(0, 0, 0))
  first <- direct_audit(z, maxitr = 0, lambda = c(.03, .1, .5), penalty = 0)
  later <- direct_audit(z, maxitr = 3, lambda = c(.03, .1, .5), penalty = 7)
  expect_equal(later$lambda.init, first$lambda.init)
  expect_length(direct_events(later, "initialize"), 1L)
  expect_equal(direct_events(later, "initialize")[[1]]$initialization.count, 1L)
  expect_equal(drop(direct_events(later, "initialize")[[1]]$grid.targets),
               drop(direct_events(first, "initialize")[[1]]$grid.targets))
  # A tiny scalar ridge remains a valid solve when its diagonal inflation
  # rounds away. It must not trigger an unnecessary 0/0 GCV calculation.
  diagonal <- list(x = matrix(c(-1, 0, 1)), B = matrix(1),
                   a = c(0, 1, 2), r = c(1, 2, 3), bw = .001)
  scalar <- direct_audit(diagonal, lambda = 1e-20, maxitr = 0)
  expect_true(all(is.finite(scalar$W)))
  expect_equal(scalar$lambda.init, 1e-20)
  expect_length(direct_events(scalar, "initialize")[[1]]$GCV, 0L)
})

test_that("direct accepted blocks satisfy independent monotone Armijo inequalities", {
  z <- direct_fixture()
  for (tau in c(1e-4, 1000)) {
    fit <- direct_audit(z, maxitr = 2, initial.w = z$w, penalty = 3,
                        control = list(tau = tau, gtol = 1e-12))
    trials <- direct_events(fit, "trial")
    expect_gt(length(trials), 0L)
    accepted <- Filter(function(item) isTRUE(item$accepted), trials)
    expect_gt(length(accepted), 0L)
    if (tau == 1e-4) expect_true(trials[[1]]$accepted)
    if (tau == 1000) expect_true(any(!vapply(trials, function(item) item$accepted, logical(1))))
    for (item in trials) {
      start <- direct_reference(item$B.start, drop(item$W.start), z$x, z$a, z$r, z$bw, 3)$fn
      direct_expect_close(item$fn.start, start)
      if (isTRUE(item$finite)) {
        trial <- direct_reference(item$B, drop(item$W), z$x, z$a, z$r, z$bw, 3)$fn
        direct_expect_close(item$fn, trial)
      }
      if (isTRUE(item$accepted)) {
        expect_true(item$finite)
        expect_lt(item$slope, 0)
        expect_lte(item$fn, start + 1e-4 * item$step * item$slope + 1e-12 * (1 + abs(start)))
        expect_lte(item$fn, start + 1e-12 * (1 + abs(start)))
        if (item$block == "W") {
          gw <- direct_reference_w(item$B.start, drop(item$W.start), z$x, z$a, z$r, z$bw, 3)
          H <- direct_reference_hessian(item$B.start, drop(item$W.start), z$x, z$a, z$r, z$bw, 3)
          expect_identical(item$direction.type, "newton")
          direction <- -drop(solve(H, gw))
          direct_expect_close(item$direction.W, direction)
          direct_expect_close(item$slope, sum(gw * direction))
          direct_expect_close(item$B, item$B.start)
          direct_expect_close(drop(item$W), drop(item$W.start) + item$step * direction)
        } else {
          G <- direct_reference_B(item$B.start, drop(item$W.start), z$x, z$a, z$r, z$bw)
          velocity <- -(G - item$B.start %*% crossprod(G, item$B.start))
          direct_expect_close(item$slope, sum(G * velocity))
          direct_expect_close(item$W, item$W.start)
        }
      }
    }
    blocks <- split(trials, vapply(trials, function(item) paste(item$sweep, item$block), character(1)))
    expect_true(all(lengths(blocks) <= 30L))
    expect_true(all(vapply(blocks, function(items)
      sum(vapply(items, function(item) item$accepted, logical(1))) <= 1L, logical(1))))
    # Each block carries its own accepted step forward. This rules out
    # secant proposals from a changing objective between the blocks.
    for (block in c("W", "B")) {
      cap <- if (block == "W") 1 else 1e10
      starts <- split(direct_events(fit, "trial", block),
                      vapply(direct_events(fit, "trial", block), function(item) item$sweep, numeric(1)))
      next.step <- min(tau, cap)
      for (items in starts) {
        direct_expect_close(items[[1]]$step, next.step)
        if (length(items) > 1L) for (k in 2:length(items))
          direct_expect_close(items[[k]]$step, .2 * items[[k - 1L]]$step)
        successes <- Filter(function(item) item$accepted, items)
        if (length(successes)) next.step <- min(successes[[1]]$step / .2, cap)
      }
    }
    expect_length(direct_events(fit, "commit"), fit$itr)
    direct_expect_returned_state(fit, z, penalty = 3, gtol = 1e-12)
  }
})

test_that("indefinite and singular coefficient curvature falls back to its true negative gradient", {
  indefinite <- direct_fixture(p = 2L)
  indefinite$r <- 100 * indefinite$r
  singular <- list(x = matrix(c(-1, -1, 1, 1)), a = c(0, 2, 0, 2),
                   r = c(0, 2, 0, 2), B = matrix(1), w = rep(.6, 4), bw = 1)
  for (z in list(indefinite, singular)) {
    H <- direct_reference_hessian(z$B, z$w, z$x, z$a, z$r, z$bw, 0)
    eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    if (nrow(z$x) == 4L) {
      expect_gte(min(eigenvalues), -1e-12)
      expect_lt(qr(H)$rank, nrow(H))
    } else {
      expect_lt(min(eigenvalues), -.1)
    }
    gw <- direct_reference_w(z$B, z$w, z$x, z$a, z$r, z$bw, 0)
    expect_gt(sqrt(sum(gw^2)), 1e-6)
    fit <- direct_audit(z, initial.w = z$w, penalty = 0,
                        control = list(tau = 1, gtol = 1e-12))
    trials <- direct_events(fit, "trial", "W")
    expect_gt(length(trials), 0L)
    expect_true(all(vapply(trials, function(item) identical(item$direction.type, "gradient"), logical(1))))
    for (item in trials) {
      direct_expect_close(item$direction.W, -gw)
      direct_expect_close(item$slope, -sum(gw^2))
      direct_expect_close(item$W, z$w - item$step * gw)
      expect_lte(item$step, 1)
      if (item$accepted) expect_lte(item$fn, item$fn.start + 1e-4 * item$step * item$slope)
    }
    expect_true(any(vapply(trials, function(item) item$accepted, logical(1))))
    if (nrow(z$x) != 4L) expect_true(any(!vapply(trials, function(item) item$accepted, logical(1))))
    direct_expect_returned_state(fit, z, penalty = 0, gtol = 1e-12)
  }
})

test_that("direction-search failure rolls back a provisional coefficient improvement atomically", {
  z <- direct_fixture()
  fit <- direct_audit(z, initial.w = z$w, B.budget = 0L,
                      control = list(tau = 1e-4, gtol = 1e-12))
  wtrials <- direct_events(fit, "trial", "W")
  expect_true(any(vapply(wtrials, function(item) item$accepted, logical(1))))
  expect_true(any(vapply(wtrials, function(item) item$accepted && item$fn < item$fn.start, logical(1))))
  expect_length(direct_events(fit, "rollback"), 1L)
  expect_length(direct_events(fit, "commit"), 0L)
  expect_equal(fit$itr, 0L)
  expect_identical(fit$stop.reason, "B_line_search_failed")
  expect_identical(fit$converge, FALSE)
  expect_equal(fit$B, z$B, tolerance = 0)
  expect_equal(drop(fit$W), z$w, tolerance = 0)
  direct_expect_returned_state(fit, z, gtol = 1e-12)
})

test_that("actual nonfinite and exhausted coefficient trials preserve the committed state", {
  z <- direct_fixture()
  w <- rep(1e-308, nrow(z$x))
  state <- direct_native(z, w = w, r = rep(-2, nrow(z$x)), penalty = 1e308)
  expect_true(state$valid)
  expect_false(state$hessian.valid)
  fit <- direct_audit(z, initial.w = w, penalty = 1e308,
                      r = rep(-2, nrow(z$x)), control = list(tau = 1e10))
  trials <- direct_events(fit, "trial", "W")
  expect_length(trials, 30L)
  expect_true(any(!vapply(trials, function(item) item$finite, logical(1))))
  expect_false(any(vapply(trials, function(item) item$accepted, logical(1))))
  expect_true(all(vapply(trials, function(item) identical(item$direction.type, "gradient"), logical(1))))
  direct_expect_close(trials[[1]]$direction.W, -drop(state$gW))
  expect_equal(trials[[1]]$step, 1)
  expect_identical(fit$stop.reason, "w_line_search_failed")
  expect_equal(fit$itr, 0L)
  expect_equal(drop(fit$W), w, tolerance = 0)
  expect_equal(fit$B, z$B, tolerance = 0)
  expect_true(all(is.finite(unlist(fit[c("fn", "value", "penalty", "gradient.norm")]))))
})

test_that("a nonfinite provisional gradient rolls back before a sweep is committed", {
  z <- list(x = matrix(c(-1, 0, 0, 1)), a = c(-1, -.01, .01, 1),
            r = c(0, 0, 1e307, 0), B = matrix(1), bw = .01)
  w <- c(-100, -50, -50, -100)
  fit <- direct_audit(z, initial.w = w, penalty = .5, control = list(tau = 1))
  accepted <- Filter(function(item) isTRUE(item$accepted), direct_events(fit, "trial", "W"))
  expect_length(accepted, 1L)
  expect_lt(accepted[[1]]$fn, accepted[[1]]$fn.start)
  expect_identical(fit$stop.reason, "numerical_failure")
  expect_false(fit$converge)
  expect_equal(fit$itr, 0L)
  expect_length(direct_events(fit, "commit"), 0L)
  expect_equal(drop(fit$W), w, tolerance = 0)
  expect_equal(fit$B, z$B, tolerance = 0)
  expect_true(all(is.finite(unlist(fit[c("fn", "value", "penalty", "gradient.norm")]))))
})

test_that("joint residuals alone certify convergence including the last budgeted sweep", {
  z <- direct_fixture()
  flat <- rep(-3, nrow(z$x))
  stationary <- direct_audit(z, r = flat, initial.w = z$w, penalty = 0, maxitr = 0)
  expect_equal(stationary$itr, 0L)
  expect_identical(stationary$stop.reason, "gradient_tolerance")
  direct_expect_returned_state(stationary, z, penalty = 0, r = flat)
  limited <- direct_audit(z, initial.w = z$w, maxitr = 0)
  expect_equal(limited$itr, 0L)
  expect_identical(limited$stop.reason, "iteration_limit")
  direct_expect_returned_state(limited, z)
  last <- direct_audit(z, r = flat, initial.w = z$w, penalty = .5,
                       maxitr = 1, control = list(tau = 1, gtol = 1e-10))
  expect_equal(last$itr, 1L)
  expect_identical(last$stop.reason, "gradient_tolerance")
  expect_lt(sqrt(sum(last$W^2)), 1e-10)
  expect_length(direct_events(last, "trial", "B"), 0L)
  direct_expect_returned_state(last, z, penalty = .5, r = flat, gtol = 1e-10)
  one <- direct_audit(z, maxitr = 1, initial.w = z$w,
                      control = list(tau = 1e-4, gtol = 1e-12))
  expect_equal(one$itr, 1L)
  expect_identical(one$stop.reason, "iteration_limit")
  direct_expect_returned_state(one, z, gtol = 1e-12)
  # Positive penalization must move beyond a nonstationary ridge initializer.
  shrink <- direct_audit(z, r = flat, penalty = .5, maxitr = 1,
                         control = list(tau = 1, gtol = 1e-10))
  expect_lt(sqrt(sum(shrink$W^2)), 1e-10)
  expect_true(shrink$converge)
})

test_that("unrepresentable steps and five small sweeps report unsuccessful stagnation", {
  z <- direct_fixture()
  tiny <- direct_audit(z, initial.w = z$w, maxitr = 8,
                       control = list(tau = 1e-300, gtol = 1e-20))
  expect_identical(tiny$stop.reason, "stagnation")
  expect_false(tiny$converge)
  expect_equal(tiny$itr, 0L)
  expect_equal(drop(tiny$W), z$w, tolerance = 0)
  direct_expect_returned_state(tiny, z, gtol = 1e-20)
  small <- direct_audit(z, initial.w = z$w, maxitr = 8,
                        control = list(tau = 1e-12, gtol = 1e-20, btol = 1e-6, ftol = 1e-6))
  expect_equal(small$itr, 5L)
  expect_identical(small$stop.reason, "stagnation")
  expect_false(small$converge)
  expect_gt(max(small$gradient.norm), 1e-6)
  direct_expect_returned_state(small, z, gtol = 1e-20)
})

test_that("legacy and optimized direct policies retain independent continuous predictions", {
  z <- direct_fixture()
  legacy <- list(B = z$B, W = matrix(rep(10, nrow(z$x))), x = z$x,
                 bw = z$bw, keep.data = TRUE, method = "direct")
  class(legacy) <- c("orthoDr", "fit", "pdose")
  fit <- direct_public(z, maxitr = 2, K = 5, keep.data = TRUE)
  one.grid <- direct_public(z, maxitr = 1, K = 1)
  direct_expect_returned_state(one.grid, z)
  for (object in list(legacy, fit)) for (n in c(1L, 3L)) {
    tx <- z$x[seq_len(n), , drop = FALSE]
    reference <- direct_reference_predict(object$B, drop(object$W), z$x, tx, z$bw)
    direct_expect_close(predict(object, tx)$pred, reference)
    restored <- unserialize(serialize(object, NULL))
    expect_identical(predict(restored, tx), predict(object, tx))
  }
  expect_gt(as.vector(predict(legacy, z$x[1, , drop = FALSE])$pred), max(z$a))
  reference <- direct_reference(fit$B, drop(fit$W), z$x, z$a, z$r, z$bw)
  direct_expect_close(predict(fit, z$x)$pred, reference$dose)
  direct_expect_returned_state(fit, z)
  no.data <- direct_public(z, maxitr = 0)
  expect_error(predict(no.data, z$x[1, , drop = FALSE]), "keep.data")
})

test_that("direct supplied controls, ridge/reward penalties and budgets are validated", {
  z <- direct_fixture()
  for (bad in list(-1, Inf, NA_real_, c(1, 2), "1"))
    expect_error(direct_public(z, lambda.reward = bad, maxitr = 0), "lambda.reward")
  for (bad in list(numeric(), 0, -.1, Inf, NA_real_, c(.1, NA_real_), "1"))
    expect_error(direct_public(z, lambda = bad, maxitr = 0), "lambda")
  for (bad in list(-1, Inf, NA_real_, c(1, 2), "1", .Machine$integer.max + 1))
    expect_error(direct_public(z, maxitr = bad), "maxitr")
  for (name in c("rho", "eta", "tau", "gtol", "btol", "ftol")) {
    for (bad in list(NULL, 0, -1, Inf, NA_real_, c(.1, .2), "bad")) {
      control <- setNames(list(bad), name)
      expect_error(direct_public(z, maxitr = 0, control = control), name)
    }
  }
  expect_error(direct_public(z, maxitr = 0, control = list(rho = 1)), "rho")
  expect_error(direct_public(z, maxitr = 0, control = list(eta = 1)), "eta")
  a <- direct_public(z, maxitr = 1.9)
  b <- direct_public(z, maxitr = 1)
  expect_equal(a, b)
})

test_that("direct rejects invalid data, bases, scales and bandwidths before fitting", {
  z <- direct_fixture()
  for (name in c("x", "a", "r")) {
    bad <- z
    bad[[name]][1] <- Inf
    expect_error(direct_public(bad, maxitr = 0), "finite|standard deviation")
    bad[[name]][1] <- NA_real_
    expect_error(direct_public(bad, maxitr = 0), "finite|standard deviation")
  }
  for (bad in list(0, -.1, Inf, NA_real_, c(.1, .2), "bad")) {
    zbad <- z
    zbad$bw <- bad
    expect_error(direct_public(zbad, maxitr = 0), "bw|bandwidth")
  }
  dependent <- z
  dependent$B[, 2] <- dependent$B[, 1]
  expect_error(direct_public(dependent, maxitr = 0), "rank|independent|orthogon|basis|B.initial")
  nonfinite <- z
  nonfinite$B[1, 1] <- Inf
  expect_error(direct_public(nonfinite, maxitr = 0), "finite|B.initial")
  zero.scale <- z
  zero.scale$x[, ] <- 1
  expect_error(direct_public(zero.scale, maxitr = 0), "scale|variance|standard deviation")
  expect_error(orthoDr_pdose(z$x, z$a, z$r, ndr = 0, B.initial = z$B, maxitr = 0), "ndr|Dimension")
  expect_error(orthoDr_pdose(z$x, z$a, z$r, ndr = ncol(z$x) + 1, maxitr = 0), "ndr")
})

test_that("an omitted direct reward penalty preserves the pseudo-direct call path", {
  z <- direct_fixture(p = 3L, d = 1L)
  pseudo <- function(...) orthoDr_pdose(z$x, z$a, z$r, ndr = 1,
    B.initial = z$B, bw = z$bw, K = 5, method = "pseudo_direct", maxitr = 1, ncore = 1, ...)
  expect_error(pseudo(lambda.reward = .1), "lambda.reward.*direct|direct.*lambda.reward")
  fit <- pseudo()
  expect_null(fit$lambda.reward)
  expect_null(fit$objective.version)
  # The original pseudo-direct ridge rule remains a final-state construction.
  reference <- direct_initial_reference(c(z[c("x", "a", "r", "bw")], list(B = fit$B)))
  expect_equal(drop(fit$W), reference$W, tolerance = 1e-9)
})

test_that("the original seed113 direct policy reaches joint stationarity and tolerates tiny start changes", {
  # Exact original tutorial data and saved initial basis. Keep the target
  # subtraction in this order so the stored investigation fixture is reproduced.
  set.seed(113)
  x <- matrix(runif(150 * 3, -1, 1), ncol = 3)
  a <- runif(150, 0, 2)
  target <- 1 + .3 * x[, 2]
  r <- 4 + x[, 1] - (a - target)^2 + rnorm(150, sd = .1)
  original.rows <- matrix(runif(5 * 3, -1, 1), ncol = 3)
  set.seed(20261027)
  independent.rows <- matrix(runif(2000 * 3, -1, 1), ncol = 3)
  newx <- rbind(original.rows, independent.rows)
  B <- matrix(c(.985882801086254923, .161846390611301710,
                 -.042905108884766875), 3, 1)
  perturbation <- B * 0
  perturbation[1, 1] <- 1e-12
  bases <- list(base = B, repeat_exact = B, sign = -B,
                 plus = B + perturbation, minus = B - perturbation)
  fits <- lapply(bases, function(basis) {
    orthoDr_pdose(x, a, r, ndr = 1, B.initial = basis,
                  bw = .38883871165870765, K = 12, lambda = .1,
                  lambda.reward = .1, maxitr = 500, ncore = 1, keep.data = TRUE)
  })
  baseline <- fits$base
  baseline.pred <- as.vector(predict(baseline, newx)$pred)
  expect_identical(fits$repeat_exact, baseline)
  for (name in names(fits)) {
    fit <- fits[[name]]
    expect_true(fit$converge, info = name)
    expect_identical(fit$stop.reason, "gradient_tolerance")
    expect_lte(fit$gradient.norm[["B"]], 1e-6)
    expect_lte(fit$gradient.norm[["W"]], 1e-6)
    expect_lte(fit$itr, 500L)
    expect_equal(crossprod(fit$B), matrix(1), tolerance = 1e-10)
    prediction <- as.vector(predict(fit, newx)$pred)
    expect_true(all(is.finite(prediction)))
    expect_lte(abs(fit$fn - baseline$fn), 1e-8 * (1 + abs(baseline$fn)))
    expect_lte(max(abs(prediction - baseline.pred)), 1e-5 * diff(range(a)))
    if (name == "sign")
      expect_lte(max(abs(prediction - baseline.pred)), 1e-8 * max(1, abs(baseline.pred)))
    # Reconstruct returned diagnostics independently, not just their labels.
    z <- list(x = x, a = a, r = r, bw = .38883871165870765)
    direct_expect_returned_state(fit, z)
  }
})

test_that("a two-direction direct policy converges consistently under signed permutations and core settings", {
  # Exact seed127 behavior fixture, including correlation and unequal units.
  # This checks the solver and representation, not the chosen penalty's quality.
  set.seed(127)
  raw <- matrix(runif(60 * 4, -1, 1), ncol = 4)
  mixing <- diag(4)
  mixing[1, 2] <- .55
  mixing[2, 3] <- -.35
  units <- c(.35, 2.3, .8, 1.6)
  x <- sweep(raw %*% mixing, 2, units, "*")
  a <- runif(60, 0, 2)
  standardized <- sweep(x, 2, units, "/")
  target <- 1 + .25 * tanh(standardized[, 1]) + .2 * sin(standardized[, 2])
  background <- 4 + .35 * tanh(x[, 1] / units[1])
  r <- background - (a - target)^2 + rnorm(60, sd = .1)
  B <- qr.Q(qr(matrix(rnorm(8), 4, 2)))
  original.rows <- sweep(matrix(runif(5 * 4, -1, 1), ncol = 4) %*% mixing, 2, units, "*")
  set.seed(20261154)
  independent.rows <- sweep(matrix(runif(2000 * 4, -1, 1), ncol = 4) %*% mixing, 2, units, "*")
  newx <- rbind(original.rows, independent.rows)
  fit <- function(basis, cores) {
    orthoDr_pdose(x, a, r, ndr = 2, B.initial = basis,
                  bw = .50540723922844411, K = 12, lambda = .1,
                  lambda.reward = .1, maxitr = 500, ncore = cores, keep.data = TRUE)
  }
  fits <- list(base = fit(B, 1),
                 signed_permutation = fit(B[, 2:1, drop = FALSE] %*% diag(c(-1, 1)), 1),
                 two_cores = fit(B, 2))
  baseline <- fits$base
  baseline.pred <- as.vector(predict(baseline, newx)$pred)
  z <- list(x = x, a = a, r = r, bw = .50540723922844411)
  for (name in names(fits)) {
    candidate <- fits[[name]]
    expect_true(candidate$converge, info = name)
    expect_identical(candidate$stop.reason, "gradient_tolerance")
    expect_lte(candidate$itr, 500L)
    direct_expect_returned_state(candidate, z)
    prediction <- as.vector(predict(candidate, newx)$pred)
    expect_true(all(is.finite(prediction)))
    tolerance <- if (name == "two_cores") 1e-7 else 1e-8
    expect_lte(abs(candidate$fn - baseline$fn), tolerance * (1 + abs(baseline$fn)))
    expect_lte(max(abs(prediction - baseline.pred)), tolerance * max(1, abs(baseline.pred)))
  }
})
