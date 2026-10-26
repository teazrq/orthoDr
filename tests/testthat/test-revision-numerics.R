# Independent arithmetic references for the targeted optimizer/dose repairs.
revision_cayley <- function(B, G, tau) {
  skew <- (G %*% t(B) - B %*% t(G)) / 2
  solve(diag(nrow(B)) + tau * skew, B - tau * skew %*% B)
}

revision_dose_data <- function() {
  set.seed(86723)
  n <- 35L
  x <- matrix(rnorm(n * 3L), n, 3L)
  a <- runif(n, 0, 3)^2 / 3
  r <- 3 + x[, 1] - (a - (1 + .3 * x[, 2]))^2
  list(x = x, a = a, r = r, B = matrix(c(.6, .8, 0), 3, 1), bw = .6)
}

revision_x_kernel <- function(B, x, bw) {
  z <- x %*% B
  z <- sweep(z, 2, apply(z, 2, sd) * bw * sqrt(2), "/")
  exp(-as.matrix(dist(z))^2)
}

revision_dose_w <- function(B, d, grid, lambda = .1) {
  k <- revision_x_kernel(B, d$x, d$bw)
  a_dist <- exp(-(outer(d$a, grid, "-") / (sd(d$a) * d$bw))^2)
  means <- vapply(seq_along(grid), function(j) {
    weights <- sweep(k, 1, a_dist[, j], "*")
    drop(crossprod(weights, d$r)) / colSums(weights)
  }, numeric(nrow(d$x)))
  target <- grid[max.col(means, ties.method = "first")]
  scores <- vapply(lambda, function(penalty) {
    smoother <- k %*% solve(k + diag(penalty, nrow(k)))
    nrow(k) * sum(((diag(nrow(k)) - smoother) %*% target)^2) /
      sum(diag(diag(nrow(k)) - smoother))^2
  }, numeric(1))
  unname(drop(solve(k + diag(lambda[which.min(scores)], nrow(k)), target)))
}

revision_dose_f <- function(B, W, d) {
  k <- revision_x_kernel(B, d$x, d$bw)
  dose <- drop(crossprod(k, W))
  weights <- k * exp(-(outer(d$a, dose, "-") / (sd(d$a) * d$bw))^2)
  -mean(drop(crossprod(weights, d$r)) / colSums(weights))
}

revision_semi_f <- function(B, d) {
  z <- cbind(d$x %*% B, d$a)
  z <- sweep(z, 2, apply(z, 2, sd) * d$bw * sqrt(2), "/")
  weights <- exp(-as.matrix(dist(z))^2)
  fitted <- drop(crossprod(weights, d$r)) / colSums(weights)
  sqrt(mean((d$r - fitted)^2))
}

test_that("stationary and flat starts return zero completed updates in both solve branches", {
  for (p in c(3L, 6L)) {
    B <- matrix(c(1, rep(0, p - 1L)), p, 1)
    A <- diag(seq(p, 1L))
    fn <- function(B) -sum(B * (A %*% B)) / 2
    gr <- function(B) -A %*% B
    fit <- ortho_optim(B, fn, gr, maxitr = 8L)
    expect_identical(fit$converge, TRUE)
    expect_equal(fit$itr, 0L)
    expect_equal(fit$B, B)
    expect_equal(fit$fn, -p / 2)
    expect_equal(as.vector(fit$fn_Seq), rep(0, 8L))

    flat <- ortho_optim(B, function(B) 1L, function(B) matrix(0L, p, 1), maxitr = 8L)
    expect_identical(flat$converge, TRUE)
    expect_equal(flat$itr, 0L)
    numerical <- ortho_optim(B, function(B) matrix(1L, 1, 1), maxitr = 8L)
    expect_identical(numerical$converge, TRUE)
    expect_equal(numerical$itr, 0L)
  }
})

test_that("finite tiny steps cannot make nonstationary movement windows converge", {
  for (p in c(3L, 6L)) {
    B <- matrix(c(.6, .8, rep(0, p - 2L)), p, 1)
    A <- diag(c(4, 2, rep(1, p - 2L)))
    fn <- function(B) -sum(B * (A %*% B)) / 2
    gr <- function(B) -A %*% B
    fit <- ortho_optim(B, fn, gr, control = list(tau = 1e-20), maxitr = 8L)
    projected <- gr(fit$B) - fit$B %*% crossprod(gr(fit$B), fit$B)
    expect_gt(norm(projected, "F"), .9)
    expect_identical(fit$converge, FALSE)
    expect_equal(fit$itr, 8L)
    expect_equal(fit$B, B)
    expect_true(all(is.finite(fit$fn_Seq)))
  }
})

test_that("small valid BB ratios and denominators follow independent Cayley updates", {
  B <- matrix(c(.6, .8, 0), 3, 1)
  A <- diag(c(4, 2, 1))
  # Each scale makes one of the valid quotient and its denominator very small.
  for (scale in c(1e-12, 1e12)) {
    fn <- function(B) -scale * sum(B * (A %*% B)) / 2
    gr <- function(B) -scale * A %*% B
    tau <- 1e-4 / scale
    next_B <- revision_cayley(B, gr(B), tau)
    old_projected <- gr(B) - B %*% crossprod(gr(B), B)
    new_projected <- gr(next_B) - next_B %*% crossprod(gr(next_B), next_B)
    S <- next_B - B
    Y <- new_projected - old_projected
    SY <- abs(sum(S * Y))
    YY <- sum(Y * Y)
    candidate <- SY / YY
    expect_true(is.finite(candidate) && candidate > 0)
    if (scale < 1) expect_lt(YY, 1e-20) else expect_lt(candidate, 1e-10)

    # Independent two-step reference retains the existing five-trial search.
    Cval <- (.85 * fn(B) + fn(next_B)) / 1.85
    reference <- NULL
    for (trial in seq_len(5L)) {
      reference <- revision_cayley(next_B, gr(next_B), candidate)
      if (fn(reference) <= Cval - candidate * 1e-4 * sum(new_projected^2) || trial == 5L) break
      candidate <- candidate * .2
    }
    fit <- ortho_optim(B, fn, gr, maxitr = 2L,
                       control = list(tau = tau, gtol = 1e-30))
    expect_equal(fit$B, reference, tolerance = 1e-9)
    expect_equal(fit$itr, 2L)
    expect_equal(drop(t(fit$B) %*% fit$B), 1, tolerance = 1e-12)
  }
})

test_that("completed iteration limits and maximization histories are accurate", {
  B <- matrix(c(.6, .8, 0), 3, 1)
  A <- diag(c(4, 2, 1))
  fn <- function(B) -sum(B * (A %*% B)) / 2
  gr <- function(B) -A %*% B
  fit <- ortho_optim(B, fn, gr, maxitr = 100L)
  at_limit <- ortho_optim(B, fn, gr, maxitr = fit$itr)
  expect_identical(at_limit$converge, TRUE)
  expect_equal(at_limit$itr, fit$itr)
  expect_equal(at_limit$B, fit$B)
  exhausted <- ortho_optim(B, fn, gr, maxitr = 1L)
  expect_identical(exhausted$converge, FALSE)
  expect_equal(exhausted$itr, 1L)

  maximum <- ortho_optim(B, function(B) -fn(B), function(B) -gr(B),
                         maximize = TRUE, maxitr = 100L)
  expect_equal(maximum$fn, -fit$fn)
  expect_equal(maximum$fn_Seq, -fit$fn_Seq)
  expect_equal(as.vector(maximum$fn_Seq)[maximum$itr], maximum$fn)
})

test_that("callbacks are checked on later evaluations and accept integer storage", {
  B <- matrix(c(.6, .8, 0), 3, 1)
  A <- diag(c(4, 2, 1))
  fn <- function(B, A) -sum(B * (A %*% B)) / 2
  gr <- function(B, A) -A %*% B
  # Expressions in ... are ordinary closure arguments.
  expect_s3_class(ortho_optim(B, fn, gr, A = diag(c(4, 2, 1)), maxitr = 1), "orthoDr")
  integer_fit <- ortho_optim(matrix(c(1, 0, 0), 3, 1), function(B) 1L,
                            function(B) matrix(c(1L, 0L, 0L), 3, 1))
  expect_equal(integer_fit$itr, 0L)

  for (bad in list(numeric(0), c(1, 2), "bad", Inf, NaN, NA_integer_)) {
    calls <- 0L
    changing_fn <- function(B) {
      calls <<- calls + 1L
      if (calls > 2L) bad else fn(B, A)
    }
    expect_error(ortho_optim(B, changing_fn, function(B) gr(B, A)), "fn must return")
    expect_gt(calls, 2L)
  }
  for (bad in list(c(1, 2, 3), matrix(1, 1, 3), matrix("bad", 3, 1),
                   matrix(Inf, 3, 1), matrix(NA_integer_, 3, 1))) {
    calls <- 0L
    changing_grad <- function(B) {
      calls <<- calls + 1L
      if (calls > 2L) bad else gr(B, A)
    }
    expect_error(ortho_optim(B, function(B) fn(B, A), changing_grad), "grad must return")
    expect_gt(calls, 2L)
  }
})

test_that("dose wrapper uses one observed scale, resolves defaults and floors K once", {
  stub <- function(...) {
    args <- list(...)
    list(B = args[[1]], fn = 0, itr = 0L, converge = TRUE,
         a_dist = args[[4]], grid = args[[5]])
  }
  local_mocked_bindings(pdose_direct_solver = stub, .package = "orthoDr")
  set.seed(1026)
  x <- matrix(rnorm(150L * 3L), 150L, 3L)
  a <- seq(.1, 3, length.out = 150L)^2
  B <- matrix(c(1, 0, 0), 3, 1)
  run <- function(a, ...) orthoDr_pdose(x, a, seq_along(a), ndr = 1,
                                       B.initial = B, bw = .6, ...)
  fit <- run(a)
  explicit <- run(a, method = "direct", K = 12L)
  expect_identical(fit$method, "direct")
  expect_equal(fit$grid, explicit$grid)
  expect_equal(ncol(fit$a_dist), 12L)
  reference <- exp(-(outer(a, as.vector(fit$grid), "-") / (sd(a) * .6))^2)
  expect_equal(fit$a_dist, reference, tolerance = 1e-13)
  expect_equal(fit$a_dist[1, 1], 1)
  expect_equal(run(7 * a)$a_dist, fit$a_dist, tolerance = 1e-13)
  expect_equal(run(a + 9)$a_dist, fit$a_dist, tolerance = 1e-13)
  expect_equal(run(a, K = 5.8)$grid, run(a, K = 5)$grid)
  one <- run(a, K = 1)
  expect_equal(as.vector(one$grid), min(a))
  expect_equal(dim(one$a_dist), c(150L, 1L))
  for (bad in list(0, .9, Inf, NA_real_, c(1, 2), "3"))
    expect_error(run(a, K = bad), "K must be")
  expect_error(run(rep(1, 150)), "positive finite standard deviation")
  expect_error(run(replace(a, 1, Inf)), "positive finite standard deviation")
  expect_error(run(a, method = "unknown"), "arg")
})

test_that("direct dose matches an independent conditional-gradient update and final state", {
  d <- revision_dose_data()
  grid <- seq(min(d$a), max(d$a), length.out = 7L)
  W <- revision_dose_w(d$B, d, grid, c(.03, .1, .5))
  f0 <- revision_dose_f(d$B, W, d)
  epsilon <- 1e-6
  G <- matrix(vapply(seq_len(nrow(d$B)), function(i) {
    perturbed <- d$B
    perturbed[i, 1] <- perturbed[i, 1] + epsilon
    (revision_dose_f(perturbed, W, d) - f0) / epsilon
  }, numeric(1)), nrow(d$B), 1)
  projected <- G - d$B %*% crossprod(G, d$B)
  tau <- .001
  reference_B <- NULL
  for (trial in seq_len(5L)) {
    reference_B <- revision_cayley(d$B, G, tau)
    reference_W <- revision_dose_w(reference_B, d, grid, c(.03, .1, .5))
    reference_f <- revision_dose_f(reference_B, reference_W, d)
    if (reference_f <= f0 - tau * 1e-4 * sum(projected^2) || trial == 5L) break
    tau <- .2 * tau
  }
  fit <- orthoDr_pdose(d$x, d$a, d$r, ndr = 1, B.initial = d$B, bw = d$bw,
                       lambda = c(.03, .1, .5), K = 7, method = "direct", maxitr = 1, ncore = 1)
  expect_equal(fit$B, reference_B, tolerance = 1e-9)
  expect_equal(as.vector(fit$W), reference_W, tolerance = 1e-8)
  expect_equal(fit$fn, reference_f, tolerance = 1e-10)
  expect_equal(fit$fn, revision_dose_f(fit$B, as.vector(fit$W), d), tolerance = 1e-12)
  expect_equal(fit$itr, 1L)
  expect_identical(fit$converge, FALSE)

  D <- matrix(c(-.8, .6, 0), 3, 1)
  plus <- d$B * cos(epsilon) + D * sin(epsilon)
  minus <- d$B * cos(epsilon) - D * sin(epsilon)
  conditional <- (revision_dose_f(plus, W, d) - revision_dose_f(minus, W, d)) / (2 * epsilon)
  profile <- (revision_dose_f(plus, revision_dose_w(plus, d, grid, c(.03, .1, .5)), d) -
                revision_dose_f(minus, revision_dose_w(minus, d, grid, c(.03, .1, .5)), d)) / (2 * epsilon)
  expect_gt(abs(conditional - profile), .001)
})

test_that("both dose methods support one grid candidate and retain their own objectives", {
  d <- revision_dose_data()
  for (method in c("direct", "pseudo_direct")) {
    for (K in c(1L, 7L)) {
      fit <- orthoDr_pdose(d$x, d$a, d$r, ndr = 1, B.initial = d$B,
                           bw = d$bw, K = K, method = method, maxitr = 1, ncore = 1)
      expect_true(all(is.finite(fit$B)))
      expect_true(all(is.finite(fit$W)))
      grid <- seq(min(d$a), max(d$a), length.out = K)
      expect_equal(as.vector(fit$W), revision_dose_w(fit$B, d, grid), tolerance = 1e-9)
      reference_f <- if (method == "direct") revision_dose_f(fit$B, as.vector(fit$W), d) else revision_semi_f(fit$B, d)
      expect_equal(fit$fn, reference_f, tolerance = 1e-12)
    }
  }
  defaults <- orthoDr_pdose(d$x, d$a, d$r, ndr = 1, B.initial = d$B,
                            bw = d$bw, K = 7, maxitr = 1, ncore = 1)
  explicit <- orthoDr_pdose(d$x, d$a, d$r, ndr = 1, B.initial = d$B,
                            bw = d$bw, K = 7, method = "direct", maxitr = 1, ncore = 1)
  expect_equal(defaults, explicit)
})

test_that("built-in flat stops finalize policies and preserve direct stopping timing", {
  d <- revision_dose_data()
  for (cores in c(1L, 2L)) {
    flat_reg <- orthoDr:::local_solver(d$B, d$x, matrix(0, nrow(d$x), 1), d$bw,
                                        1e-4, .2, .85, .001, 1e-6, 1e-6, 1e-6, 1e-6, 2L, 0L, cores)
    expect_equal(flat_reg$itr, 0L)
    expect_identical(flat_reg$converge, TRUE)
    semi <- orthoDr_pdose(d$x, d$a, rep(0, nrow(d$x)), ndr = 1, B.initial = d$B,
                          bw = d$bw, K = 1, method = "pseudo_direct", maxitr = 8, ncore = cores)
    expect_equal(semi$itr, 0L)
    expect_identical(semi$converge, TRUE)
    flat_d <- d
    flat_d$r[] <- 0
    expect_equal(as.vector(semi$W), revision_dose_w(semi$B, flat_d, min(d$a)), tolerance = 1e-10)
    direct <- orthoDr_pdose(d$x, d$a, rep(0, nrow(d$x)), ndr = 1, B.initial = d$B,
                            bw = d$bw, K = 1, method = "direct", maxitr = 5, ncore = cores)
    expect_equal(direct$itr, 5L)
    expect_identical(direct$converge, TRUE)
    expect_equal(direct$B, d$B)
  }
})

test_that("built-in nonstationary stagnation respects effective limits and core agreement", {
  d <- revision_dose_data()
  for (cores in c(1L, 2L)) {
    reg <- orthoDr_reg(d$x, d$r, ndr = 1, B.initial = d$B, bw = d$bw,
                       method = "sir", maxitr = 2, ncore = cores,
                       control = list(tau = 1e-20, gtol = 1e-20))
    expect_equal(reg$itr, 5L)
    expect_identical(reg$converge, FALSE)
    for (method in c("direct", "pseudo_direct")) {
      fit <- orthoDr_pdose(d$x, d$a, d$r, ndr = 1, B.initial = d$B, bw = d$bw,
                           K = 7, method = method, maxitr = 8, ncore = cores,
                           control = list(tau = 1e-20, gtol = 1e-20))
      expect_equal(fit$itr, 8L)
      expect_identical(fit$converge, FALSE)
      expect_equal(fit$B, d$B)
    }
  }
  for (method in c("direct", "pseudo_direct")) {
    fits <- lapply(c(1L, 2L), function(cores) orthoDr_pdose(d$x, d$a, d$r,
                       ndr = 1, B.initial = d$B, bw = d$bw, K = 7,
                       method = method, maxitr = 3, ncore = cores))
    expect_equal(fits[[1]], fits[[2]], tolerance = 1e-10)
  }
})
