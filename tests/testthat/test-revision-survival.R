# Independent complete-risk-set references. Native values are obtained from
# zero-update solves; a one-update comparison also checks the native gradient.
surv_revision_preprocess <- function(x, time, event) {
  oo <- order(time)
  x <- x[oo, , drop = FALSE]
  time <- time[oo]
  event <- event[oo]
  ranks <- rank(time) / (length(time) + 1)
  yt <- if (sd(ranks) == 0) rep(0, length(time)) else
    drop(scale(ranks)) / sqrt(2) / silverman(1, length(time))
  fail <- which(event == 1)
  phi <- matrix(0, ncol(x), length(fail))
  for (h in seq_along(fail)) {
    j <- fail[h]
    weights <- event * exp(-abs(yt - yt[j])) * (time >= time[j])
    phi[, h] <- colSums(x * weights) / sum(weights) -
      colMeans(x[time >= time[j], , drop = FALSE])
  }
  list(x = x, time = time, event = event, fail = fail,
       risk = match(time[fail], time), phi = phi)
}

surv_revision_objective <- function(B, dat, method, bw = .6, legacy = FALSE) {
  x <- dat$x
  z <- x %*% B
  z <- sweep(z, 2, apply(z, 2, sd) * bw * sqrt(2), "/")
  kernel <- exp(-as.matrix(dist(z))^2)
  ee <- if (method == "forward") numeric(ncol(x)) else
    matrix(0, ncol(x), ncol(x))
  for (h in seq_along(dat$fail)) {
    j <- dat$fail[h]
    risk <- if (legacy) seq.int(j, nrow(x)) else
      which(dat$time >= dat$time[j])
    if (method == "dm") {
      for (i in risk) {
        weights <- kernel[i, risk]
        residual <- x[i, ] - colSums(x[risk, , drop = FALSE] * weights) / sum(weights)
        ee <- ee + ((i == j) - kernel[i, j] / sum(weights)) *
          outer(dat$phi[, h], residual)
      }
    } else {
      weights <- kernel[risk, j]
      residual <- x[j, ] - colSums(x[risk, , drop = FALSE] * weights) / sum(weights)
      ee <- ee + if (method == "forward") residual else outer(dat$phi[, h], residual)
    }
  }
  sum(ee^2) / length(dat$fail)^2
}

surv_revision_grouped_dm <- function(B, dat, bw = .6) {
  z <- dat$x %*% B
  z <- sweep(z, 2, apply(z, 2, sd) * bw * sqrt(2), "/")
  kernel <- exp(-as.matrix(dist(z))^2)
  ee <- matrix(0, ncol(dat$x), ncol(dat$x))
  for (time in unique(dat$time[dat$fail])) {
    failures <- which(dat$time == time & dat$event == 1)
    risk <- which(dat$time >= time)
    phi <- dat$phi[, match(failures[1], dat$fail)]
    for (i in risk) {
      weights <- kernel[i, risk]
      residual <- dat$x[i, ] -
        colSums(dat$x[risk, , drop = FALSE] * weights) / sum(weights)
      increment <- as.numeric(i %in% failures) - sum(kernel[i, failures]) / sum(weights)
      ee <- ee + increment * outer(phi, residual)
    }
  }
  sum(ee^2) / length(dat$fail)^2
}

surv_revision_native <- function(B, dat, method, maxitr = 0L,
                                 tau = .02, gtol = 1e-12, ncore = 1L,
                                 btol = 1e-12, ftol = 1e-12) {
  solver <- get(paste0("surv_", method, "_solver"), asNamespace("orthoDr"))
  args <- list(B = B, X = dat$x, Fail_Ind = dat$fail, Risk_Start = dat$risk,
               bw = .6, rho = 1e-4, eta = .2, gamma = .85,
               tau = tau, epsilon = 1e-6, btol = btol, ftol = ftol,
               gtol = gtol, maxitr = maxitr, verbose = 0L, ncore = ncore)
  if (method != "forward") args$Phit <- dat$phi
  do.call(solver, args)
}

surv_revision_step <- function(B, dat, method, tau = .02) {
  objective <- function(b) surv_revision_objective(b, dat, method)
  f0 <- objective(B)
  gradient <- B * 0
  for (j in seq_len(ncol(B))) for (i in seq_len(nrow(B))) {
    trial <- B
    trial[i, j] <- trial[i, j] + 1e-6
    gradient[i, j] <- (objective(trial) - f0) / 1e-6
  }
  tangent <- gradient - B %*% crossprod(gradient, B)
  for (ls in seq_len(5)) {
    if (ncol(B) < floor(nrow(B) / 2)) {
      u <- cbind(gradient, B)
      v <- cbind(B, -gradient)
      trial <- B - u %*% (tau * solve(diag(2 * ncol(B)) +
        .5 * tau * crossprod(v, u), crossprod(v, B)))
    } else {
      h <- .5 * (tcrossprod(gradient, B) - tcrossprod(B, gradient))
      trial <- solve(diag(nrow(B)) + tau * h, B - tau * h %*% B)
    }
    if (objective(trial) <= f0 - tau * 1e-4 * sum(tangent^2) || ls == 5)
      return(trial)
    tau <- .2 * tau
  }
}

test_that("survival native objectives and gradient updates use complete risk sets", {
  set.seed(881)
  for (p in c(3L, 6L)) { # Both feasible-update linear-solve branches.
    x <- matrix(rnorm(24 * p), 24, p)
    time <- rep(seq_len(8), each = 3)
    event <- rep(c(1, 0, 1), 8)
    B <- qr.Q(qr(matrix(rnorm(p), p, 1)))
    perm <- as.vector(rbind(seq(3, 24, 3), seq(1, 24, 3), seq(2, 24, 3)))
    dat <- surv_revision_preprocess(x, time, event)
    shuffled <- surv_revision_preprocess(x[perm, , drop = FALSE], time[perm], event[perm])
    for (method in c("forward", "dn", "dm")) {
      expected <- surv_revision_objective(B, dat, method)
      if (method == "dm")
        expect_equal(expected, surv_revision_grouped_dm(B, dat), tolerance = 1e-11)
      fit <- surv_revision_native(B, dat, method)
      expect_equal(fit$fn, expected, tolerance = 1e-11)
      expect_equal(surv_revision_native(B, shuffled, method)$fn, expected, tolerance = 1e-11)
      nextB <- surv_revision_native(B, dat, method, maxitr = 1L)$B
      expect_equal(nextB, surv_revision_step(B, dat, method), tolerance = 2e-7)
      expect_equal(surv_revision_native(B, shuffled, method, maxitr = 1L)$B,
                   nextB, tolerance = 2e-7)
      expect_equal(surv_revision_native(B, dat, method, ncore = 2L)$fn,
                   expected, tolerance = 1e-11)
    }
    untied <- surv_revision_preprocess(x, seq_len(24), event)
    tied <- surv_revision_preprocess(x, rep(1, 24), event)
    for (method in c("forward", "dn", "dm")) {
      expect_equal(surv_revision_native(B, untied, method)$fn,
                   surv_revision_objective(B, untied, method, legacy = TRUE), tolerance = 1e-11)
      expect_equal(surv_revision_native(B, tied, method)$fn,
                   surv_revision_objective(B, tied, method), tolerance = 1e-11)
      if (method == "dm")
        expect_equal(surv_revision_native(B, tied, method)$fn,
                     surv_revision_grouped_dm(B, tied), tolerance = 1e-11)
    }
  }
})

surv_revision_cp_moment <- function(x, time, event, bw) {
  p <- ncol(x)
  eig <- eigen(cov(x))
  inverse <- eig$vectors %*% diag(1 / sqrt(eig$values), p, p) %*% t(eig$vectors)
  white <- scale(x, scale = FALSE) %*% inverse
  fail <- which(event == 1)[order(time[event == 1])]
  times <- time[fail]
  moment <- matrix(0, p, p)
  for (j in seq_along(fail)) {
    midpoint <- mean(which(times == times[j]))
    lower <- floor(max(1, midpoint - length(fail) * bw / 2))
    upper <- floor(min(length(fail), midpoint + length(fail) * bw / 2))
    neighborhood <- event == 1 & time >= times[lower] & time <= times[upper]
    riskmean <- colMeans(white[time >= times[j], , drop = FALSE])
    moment <- moment + outer(white[fail[j], ] - riskmean,
                             colMeans(white[neighborhood, , drop = FALSE]) - riskmean)
  }
  list(moment = moment, inverse = inverse)
}

test_that("CP-SIR tied neighborhoods match an independent estimation matrix", {
  set.seed(411)
  x <- matrix(rnorm(30 * 3), 30, 3)
  time <- rep(seq_len(10), each = 3)
  event <- rep(c(1, 1, 0), 10)
  perm <- as.vector(rbind(seq(3, 30, 3), seq(1, 30, 3), seq(2, 30, 3)))
  for (times in list(time, seq_len(30), rep(1, 30))) {
    direct <- surv_revision_cp_moment(x, times, event, bw = .35)
    fit <- CP_SIR(x, times, event, bw = .35)
    other <- CP_SIR(x[perm, , drop = FALSE], times[perm], event[perm], bw = .35)
    decomposition <- svd(direct$moment)
    expect_equal(fit$values, decomposition$d, tolerance = 1e-10)
    expect_equal(other$values, fit$values, tolerance = 1e-10)
    expect_gt(decomposition$d[1] - decomposition$d[2], 1e-3)
    leading <- direct$inverse %*% decomposition$v[, 1, drop = FALSE]
    leading <- leading / sqrt(sum(leading^2))
    expect_equal(tcrossprod(fit$vectors[, 1]), tcrossprod(leading[, 1]), tolerance = 1e-9)
    expect_equal(tcrossprod(other$vectors[, 1]), tcrossprod(fit$vectors[, 1]), tolerance = 1e-9)
  }
})

test_that("CP-SIR preserves one-feature and one-event matrices without a moment-rank gate", {
  x <- matrix(2 * seq_len(12), 12, 1)
  expect_equal(dim(CP_SIR(x, seq_len(12), rep(c(1, 0), 6))$vectors), c(1L, 1L))
  set.seed(718)
  x <- matrix(rnorm(18 * 3), 18, 3)
  for (index in c(1L, 18L)) {
    event <- as.numeric(seq_len(18) == index)
    cp <- CP_SIR(x, seq_len(18), event)
    expect_equal(dim(cp$vectors), c(3L, 3L))
    expect_true(all(is.finite(cp$vectors)))
    fit <- orthoDr_surv(x, seq_len(18), event, method = "dn", ndr = 2,
                       maxitr = 0, ncore = 1)
    expect_equal(dim(fit$B), c(3L, 2L))
    expect_true(is.finite(fit$fn))
  }
  expect_error(CP_SIR(x, seq_len(18), rep(0, 18)), "at least one event")
  expect_error(orthoDr_surv(x, seq_len(18), rep(0, 18)), "at least one event")
  singular <- cbind(seq_len(18), seq_len(18))
  event <- rep(c(1, 0), 9)
  expect_error(CP_SIR(singular, seq_len(18), event), "positive definite")
  expect_no_error(orthoDr_surv(singular, seq_len(18), event, method = "dn", ndr = 1,
                             B.initial = matrix(c(1, 0), 2, 1), maxitr = 0, ncore = 1))
})

test_that("public survival fits retain ties, dimensions, and bandwidth for every method", {
  set.seed(842)
  x <- matrix(rnorm(24 * 3), 24, 3)
  times <- rep(seq_len(8), each = 3)
  event <- rep(c(1, 0, 1), 8)
  perm <- as.vector(rbind(seq(3, 24, 3), seq(1, 24, 3), seq(2, 24, 3)))
  for (method in c("forward", "dn", "dm")) {
    for (time in list(times, rep(1, 24))) for (initialize in c("default", "explicit")) {
      basis <- if (initialize == "default") NULL else matrix(c(1, 0, 0), 3, 1)
      fit <- orthoDr_surv(x, time, event, method = method, ndr = 1, B.initial = basis,
                         bw = .6, keep.data = TRUE, maxitr = 1, ncore = 1)
      other <- orthoDr_surv(x[perm, , drop = FALSE], time[perm], event[perm],
                           method = method, ndr = 1, B.initial = basis,
                           bw = .6, maxitr = 1, ncore = 1)
      expect_equal(fit$bw, .6)
      expect_equal(dim(fit$B), c(3L, 1L))
      expect_true(is.finite(fit$fn))
      expect_equal(other$fn, fit$fn, tolerance = 1e-7)
      expect_equal(tcrossprod(other$B), tcrossprod(fit$B), tolerance = 1e-7)
    }
    one <- orthoDr_surv(matrix(3 * seq_len(12), 12, 1), seq_len(12),
                        rep(c(1, 0), 6), method = method, ndr = 1, maxitr = 1, ncore = 1)
    expect_equal(dim(one$B), c(1L, 1L))
    expect_true(is.finite(one$fn))
    expect_equal(one$itr, 0L)
    expect_true(one$converge)
  }
})

test_that("survival loops report completed updates and reject false stagnation success", {
  set.seed(413)
  x <- matrix(rnorm(24 * 3), 24, 3)
  dat <- surv_revision_preprocess(x, rep(seq_len(8), each = 3), rep(c(1, 0, 1), 8))
  B <- matrix(c(1, 2, 3) / sqrt(14), 3, 1)
  for (method in c("forward", "dn", "dm")) {
    no_updates <- surv_revision_native(B, dat, method)
    expect_equal(no_updates$itr, 0L)
    expect_false(no_updates$converge)
    exhausted <- surv_revision_native(B, dat, method, maxitr = 1L)
    expect_equal(exhausted$itr, 1L)
    expect_false(exhausted$converge)
    stagnant <- surv_revision_native(B, dat, method, maxitr = 7L, tau = 1e-300,
                                     btol = 1e-8, ftol = 1e-8)
    expect_false(stagnant$converge)
    expect_equal(stagnant$itr, 7L)
    expect_true(is.finite(stagnant$fn))
  }
})
