spectral_fixture <- function(p = 3L, d = 2L) {
  n <- 12L
  index <- seq_len(n * p)
  x <- matrix(sin(sqrt(2) * index) + cos(sqrt(3) * index) + .1 * sin(index^2 / 7), n, p)
  B <- qr.Q(qr(matrix(sin(seq_len(p * d)), p, d)))
  a <- seq(.1, 2.3, length.out = n)
  r <- 3 + .2 * x[, 1] - (a - 1.1 - .3 * tanh(x[, min(2L, p)]))^2
  list(x = x, a = a, r = r, B = B, w = .13 + .09 * cos(seq_len(n)),
       alpha = 1.1, bw = .8)
}

spectral_fit <- function(z, maxitr = 6L, intercept = TRUE, fail = -1L,
                         tau = .001, gtol = 1e-12, penalty = .17) {
  grid <- seq(min(z$a), max(z$a), length.out = 7L)
  orthoDr:::pdose_direct_audit_solver(z$B, z$x, z$a,
    exp(-outer(z$a, grid, "-")^2 / (sd(z$a) * z$bw)^2), grid,
    z$r, .1, z$bw, 1e-4, .2, .85, tau, 1e-6, 1e-10, 1e-10,
    gtol, as.integer(maxitr), 0L, 1L, penalty, TRUE, 30L, z$w,
    intercept, if (intercept) range(z$a) else NULL,
    if (intercept) z$alpha else NULL, as.integer(fail))
}

spectral_events <- function(fit, event, block = NULL) {
  Filter(function(e) identical(e$event, event) &&
           (is.null(block) || identical(e$block, block)), fit$trace)
}

spectral_project <- function(B, Z) {
  Z - B %*% ((crossprod(B, Z) + crossprod(Z, B)) / 2)
}

# Independent equations for the basis gradient, including the moving kernel
# means. The secant tests never call the native objective/gradient evaluator.
spectral_reference_gradient <- function(z, B, w, alpha) {
  w <- drop(w)
  n <- nrow(z$x)
  u <- z$x %*% B
  variance <- vapply(seq_len(ncol(B)), function(k) var(u[, k]), numeric(1))
  h2 <- var(z$a) * z$bw^2
  logK <- matrix(0, n, n)
  for (k in seq_len(ncol(B)))
    logK <- logK - outer(u[, k], u[, k], "-")^2 / (2 * z$bw^2 * variance[k])
  K <- exp(logK)
  dose <- drop(sweep(K, 2, colMeans(K), "-") %*% w) + alpha
  logits <- logK + outer(2 * z$a / h2, dose) - z$a^2 / h2
  P <- exp(sweep(logits, 2, apply(logits, 2, max), "-"))
  P <- sweep(P, 2, colSums(P), "/")
  reward <- drop(crossprod(P, z$r))
  influence <- P * sweep(matrix(z$r, n, n), 2, reward, "-")
  gdose <- -2 * colSums(influence * z$a) / (n * h2)
  sensitivity <- -influence / n + K * outer(w, gdose - mean(gdose))
  covariance <- cov(z$x)
  gradient <- matrix(0, nrow(B), ncol(B))
  for (k in seq_len(ncol(B))) {
    delta <- outer(u[, k], u[, k], "-")
    cb <- drop(covariance %*% B[, k])
    for (j in seq_len(nrow(B))) {
      derivative <- -delta * outer(z$x[, j], z$x[, j], "-") / (z$bw^2 * variance[k]) +
        delta^2 * cb[j] / (z$bw^2 * variance[k]^2)
      gradient[j, k] <- sum(sensitivity * derivative)
    }
  }
  gradient
}

test_that("spectral scalars preserve parity, signed curvature and exact fallback", {
  threshold <- (sqrt(.Machine$double.eps) * sqrt(4)) * sqrt(9)
  cases <- list(
    list("tau1 odd", 4, 2, 3, 1L, .003, TRUE, TRUE, TRUE, 2),
    list("tau2 even", 4, 2, 3, 2L, .003, TRUE, TRUE, TRUE, 2/3),
    list("negative curvature", 4, -2, 3, 1L, .003, TRUE, TRUE, FALSE, .003),
    list("zero curvature", 4, 0, 3, 1L, .003, TRUE, TRUE, FALSE, .003),
    list("zero displacement", 0, 1, 3, 1L, .003, TRUE, TRUE, FALSE, .003),
    list("zero gradient difference", 4, 1, 0, 1L, .003, TRUE, TRUE, FALSE, .003),
    list("strict threshold", 4, threshold, 9, 1L, .003, TRUE, TRUE, FALSE, .003),
    list("invalid vectors", 4, 2, 3, 1L, 1e-250, FALSE, TRUE, FALSE, 1e-250),
    list("invalid evaluation", 4, 2, 3, 1L, 1e-250, TRUE, FALSE, FALSE, 1e-250),
    list("nonfinite pair", Inf, 2, 3, 1L, .003, TRUE, TRUE, FALSE, .003),
    list("tau1 has no lower floor", 1e-60, 1e-30, 1, 1L, .003, TRUE, TRUE, TRUE, 1e-30),
    list("tau2 has no lower floor", 1, 1e30, 1e60, 2L, .003, TRUE, TRUE, TRUE, 1e-30),
    list("upper cap", 1e30, 1e15, 1, 1L, .003, TRUE, TRUE, TRUE, 1e10),
    list("overflow safe threshold", 1e308, 1e308, 1e308, 1L, .003, TRUE, TRUE, TRUE, 1),
    list("overflowing ratio", 1e308, 1e-6, 1e-320, 1L, .003, TRUE, TRUE, FALSE, .003))
  for (case in cases) {
    out <- do.call(orthoDr:::pdose_direct_spectral_scalar,
      setNames(case[2:8], c("SS", "SZ", "ZZ", "m", "fallback", "finite_vectors", "evaluation_valid")))
    expect_identical(out$usable, case[[9]], info = case[[1]])
    expect_true(is.finite(out$step), info = case[[1]])
    expect_lte(abs(out$step - case[[10]]) / max(abs(case[[10]]), .Machine$double.xmin),
               1e-14, label = case[[1]])
    if (!case[[9]]) expect_identical(out$step, case[[6]], info = case[[1]])
  }
})

test_that("same-policy tangent secants use committed pre-move anchors", {
  z <- spectral_fixture()
  fit <- spectral_fit(z)
  proposals <- spectral_events(fit, "B.proposal")
  expect_gte(length(proposals), 3L)
  expect_identical(proposals[[1]]$origin, "fallback")
  expect_false(proposals[[1]]$has.anchor)
  expect_equal(proposals[[1]]$m, 0)
  m <- 0L; anchor <- NULL; anchor.w <- NULL; anchor.alpha <- NULL
  fallback <- .001
  wrong.policy <- raw.transport <- 0
  for (entry in fit$trace) {
    if (entry$event == "B.proposal") {
      expect_equal(entry$m, m)
      expect_identical(entry$has.anchor, !is.null(anchor))
      expect_equal(entry$fallback, fallback, tolerance = 0)
      if (entry$attempted.pair) {
        expect_equal(entry$anchor, anchor, tolerance = 0)
        now <- spectral_project(entry$B.current, spectral_reference_gradient(
          z, entry$B.current, entry$W.current, entry$intercept.current))
        previous.g <- spectral_reference_gradient(z, anchor, entry$W.current, entry$intercept.current)
        previous <- spectral_project(anchor, previous.g)
        transported <- spectral_project(entry$B.current, previous)
        s <- spectral_project(entry$B.current, entry$B.current - anchor)
        delta <- now - transported
        SS <- sum(s^2); SZ <- sum(s * delta); ZZ <- sum(delta^2)
        threshold <- (sqrt(.Machine$double.eps) * sqrt(SS)) * sqrt(ZZ)
        usable <- SS > 0 && ZZ > 0 && SZ > threshold
        ratio <- if (m %% 2L == 1L) SS / SZ else SZ / ZZ
        usable <- usable && is.finite(ratio) && ratio > 0
        expect_equal(entry$gradient.current.tangent, now, tolerance = 1e-10)
        expect_equal(entry$gradient.previous.tangent, previous, tolerance = 1e-10)
        expect_equal(entry$gradient.previous.transported, transported, tolerance = 1e-10)
        expect_equal(entry$s, s, tolerance = 1e-12)
        expect_equal(entry$z, delta, tolerance = 1e-10)
        expect_equal(c(entry$SS, entry$SZ, entry$ZZ, entry$threshold),
                     c(SS, SZ, ZZ, threshold), tolerance = 1e-10)
        expect_identical(entry$origin == "spectral", usable)
        expect_equal(entry$step, if (usable) min(ratio, 1e10) else fallback, tolerance = 1e-9)
        stale <- spectral_reference_gradient(z, anchor, anchor.w, anchor.alpha)
        wrong.policy <- max(wrong.policy,
          norm(now - spectral_project(entry$B.current, spectral_project(anchor, stale)) - delta, "F"))
        raw.transport <- max(raw.transport,
          norm(now - spectral_project(entry$B.current, previous.g) - delta, "F"))
      }
    }
    if (entry$event == "trial" && entry$block == "B" && entry$accepted)
      fallback <- if (entry$step >= 1e10 * .2) 1e10 else entry$step / .2
    if (entry$event == "commit" && any(entry$B != entry$B.start)) {
      m <- m + 1L; anchor <- entry$B.start
      anchor.w <- entry$W.start; anchor.alpha <- entry$intercept.start
      expect_equal(entry$history.anchor, anchor, tolerance = 0)
      expect_equal(entry$history.m, m)
    }
  }
  # Both incorrect constructions must be distinguishable on this trajectory.
  expect_gt(wrong.policy, 1e-7)
  expect_gt(raw.transport, 1e-7)
  expect_equal(fit$audit.history$m, m)
  expect_equal(fit$audit.history$anchor, anchor, tolerance = 0)
  trials <- spectral_events(fit, "trial")
  groups <- split(trials, vapply(trials, function(e) paste(e$sweep, e$block), ""))
  expect_true(all(lengths(groups) <= 30L))
  for (entry in Filter(function(e) isTRUE(e$accepted), trials)) {
    expect_lt(entry$slope, 0)
    expect_lte(entry$fn, entry$fn.start + 1e-4 * entry$step * entry$slope)
  }
})

test_that("a final-gradient failure rolls back established state and spectral history", {
  z <- spectral_fixture()
  short <- spectral_fit(z, maxitr = 2L)
  rollback <- spectral_fit(z, maxitr = 3L, fail = 3L)
  expect_equal(short$itr, 2)
  expect_equal(rollback$itr, 2)
  expect_identical(rollback$stop.reason, "numerical_failure")
  expect_false(rollback$converge)
  fields <- c("B", "W", "intercept", "fn", "value", "penalty", "gradient.norm",
              "kernel.mean", "gradient.intercept", "audit.history")
  expect_identical(rollback[fields], short[fields])
  staged <- Filter(function(e) e$sweep == 3L && e$accepted,
                   spectral_events(rollback, "trial", "B"))
  expect_length(staged, 1L)
  expect_length(spectral_events(rollback, "rollback"), 1L)
  final <- spectral_events(rollback, "rollback")[[1]]
  expect_identical(final$history.anchor, short$audit.history$anchor)
  expect_identical(final$history.m, short$audit.history$m)
  z$B <- short$B; z$w <- drop(short$W); z$alpha <- short$intercept
  restart <- spectral_fit(z, maxitr = 1L)
  first <- spectral_events(restart, "B.proposal")[[1]]
  expect_equal(first$m, 0)
  expect_false(first$has.anchor)
  expect_false(first$attempted.pair)
  expect_identical(first$origin, "fallback")
})

test_that("legacy fits and skipped basis blocks accumulate no spectral history", {
  for (p in c(2L, 6L)) {
    fit <- spectral_fit(spectral_fixture(p), maxitr = 3L, intercept = FALSE)
    expect_equal(fit$itr, 3)
    expect_false(fit$audit.history$has.anchor)
    expect_equal(fit$audit.history$m, 0)
    proposals <- spectral_events(fit, "B.proposal")
    expect_gt(length(proposals), 0L)
    expect_true(all(vapply(proposals, function(e) !e$attempted.pair && e$origin == "fallback", logical(1))))
  }
  one <- spectral_fixture(1L, 1L)
  one$r[] <- -3
  fit <- spectral_fit(one, maxitr = 1L, tau = 1, gtol = 1e-10, penalty = .5)
  expect_true(fit$converge)
  expect_equal(fit$itr, 1)
  expect_length(spectral_events(fit, "skip", "B"), 1L)
  expect_false(fit$audit.history$has.anchor)
  expect_equal(fit$audit.history$m, 0)
})
