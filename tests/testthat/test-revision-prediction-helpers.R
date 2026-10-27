# Independent squared-distance calculation, without feature standardization.
revision_kernel_reference <- function(x, testx, B, bw) {
  z <- x %*% B
  ztest <- testx %*% B
  denominator <- apply(z, 2, sd) * bw * sqrt(2)
  out <- matrix(0, nrow(testx), nrow(x))
  for (i in seq_len(nrow(testx))) {
    for (j in seq_len(nrow(x))) {
      out[i, j] <- exp(-sum(((ztest[i, ] - z[j, ]) / denominator)^2))
    }
  }
  out
}

revision_fit_object <- function(x, y, B, bw, censor = NULL) {
  ans <- list(x = x, y = y, B = B, bw = bw, keep.data = TRUE)
  if (!is.null(censor)) ans$censor <- censor
  class(ans) <- c("orthoDr", "fit", if (is.null(censor)) "reg" else "surv")
  ans
}

test_that("prediction kernels use the fitted raw-coordinate projections", {
  set.seed(7101)
  x <- sweep(sweep(matrix(rnorm(80 * 4), 80, 4), 2, c(1, 9, .3, 4), "*"),
             2, c(7, -20, 5, 11), "+")
  testx <- x[c(3, 11, 27), , drop = FALSE] + .2
  B <- cbind(c(.6, .8, 0, 0), c(.2, -.1, .9, .3))
  y <- x[, 1] - x[, 2] / 3

  for (d in 1:2) {
    basis <- B[, seq_len(d), drop = FALSE]
    object <- revision_fit_object(x, y, basis, .45)
    for (nr in c(1L, 3L)) {
      tx <- testx[seq_len(nr), , drop = FALSE]
      k <- revision_kernel_reference(x, tx, basis, object$bw)
      expect_equal(orthoDr:::prediction_kernel(object, tx), k, tolerance = 1e-12)
      expect_equal(predict(object, tx)$pred,
                   drop(k %*% y / rowSums(k)), tolerance = 1e-12)
    }

    # Serialization preserves the existing fitted-object field convention.
    restored <- unserialize(serialize(object, NULL))
    expect_identical(predict(restored, testx), predict(object, testx))
  }

  unitx <- scale(x)
  unit_test <- scale(testx, center = colMeans(x), scale = apply(x, 2, sd))
  object <- revision_fit_object(unitx, y, B[, 1, drop = FALSE], .45)
  k <- revision_kernel_reference(unitx, unit_test, object$B, object$bw)
  expect_equal(predict(object, unit_test)$pred,
               drop(k %*% y / rowSums(k)), tolerance = 1e-12)
})

test_that("survival raw-coordinate bases retain the fitted kernel", {
  set.seed(7102)
  x <- sweep(sweep(matrix(rnorm(50 * 3), 50, 3), 2, c(.4, 8, 2), "*"),
             2, c(3, -7, 15), "+")
  testx <- x[1:2, , drop = FALSE] + .1
  std_basis <- qr.Q(qr(matrix(c(1, 2, -1, 2, -1, 1), 3, 2)))
  raw_basis <- sweep(std_basis, 1, apply(x, 2, sd), "/")
  raw_basis <- sweep(raw_basis, 2, sqrt(colSums(raw_basis^2)), "/")
  std_test <- scale(testx, center = colMeans(x), scale = apply(x, 2, sd))
  for (d in 1:2) {
    basis <- raw_basis[, seq_len(d), drop = FALSE]
    object <- revision_fit_object(x, seq_len(nrow(x)), basis, .6, rep(1, nrow(x)))
    fitted_kernel <- revision_kernel_reference(scale(x), std_test,
                                               std_basis[, seq_len(d), drop = FALSE], .6)
    expect_equal(orthoDr:::prediction_kernel(object, testx), fitted_kernel,
                 tolerance = 1e-12)
    expect_identical(predict(unserialize(serialize(object, NULL)), testx),
                     predict(object, testx))
  }
})

test_that("prediction diagnoses missing bandwidth and invalid projection scale", {
  x <- matrix(seq_len(10), 10, 1)
  object <- revision_fit_object(x, seq_len(10), matrix(1), .4, rep(1, 10))
  object$method <- "forward"
  object$bw <- NULL
  expect_error(predict(object, x[1, , drop = FALSE]), "original bandwidth.*refit")
  object$bw <- .4
  object$B[] <- 0
  expect_error(predict(object, x[1, , drop = FALSE]), "projection scales")
  object$B[] <- Inf
  expect_error(predict(object, x[1, , drop = FALSE]), "projection scales")
})

test_that("weighted survival agrees with independent survfit references", {
  cases <- list(
    list(time = 1:6, event = rep(1, 6), x = rep(c(-1, 1), 3)),
    list(time = c(1, 1, 2, 2, 3, 4), event = c(1, 0, 1, 1, 0, 1), x = c(-1, 2, 0, 1, 3, -.5)),
    list(time = 1:6, event = c(0, 0, 0, 0, 0, 1), x = seq(-1, 2, length.out = 6)),
    list(time = c(1, 1, 2, 2, 3, 4), event = c(1, 1, 1, 0, 0, 1), x = c(-2, -.5, 0, 2, 1, 3)),
    list(time = rep(2, 6), event = c(1, 1, 0, 1, 0, 1), x = rep(c(-1, 1), 3))
  )

  for (case in cases) {
    x <- matrix(case$x + 5, ncol = 1)
    object <- revision_fit_object(x, case$time, matrix(1), .7, case$event)
    testx <- matrix(c(5, 5.5), ncol = 1)
    weights <- revision_kernel_reference(x, testx, object$B, object$bw)
    ans <- predict(object, testx)
    event_times <- sort(case$time[case$event == 1])
    expect_identical(ans$timepoints, event_times)
    expect_equal(dim(ans$surv), c(sum(case$event), nrow(testx)))
    for (j in seq_len(nrow(testx))) {
      reference <- survival::survfit(
        survival::Surv(case$time, case$event) ~ 1, weights = weights[j, ]
      )
      reference <- summary(reference, times = event_times, extend = TRUE)$surv
      expect_equal(ans$surv[, j], unname(reference), tolerance = 1e-12)
      expect_true(all(diff(ans$surv[, j]) <= 0))
    }
    expect_equal(predict(object, testx[1, , drop = FALSE])$surv,
                 ans$surv[, 1, drop = FALSE])
    object$censor <- as.logical(case$event)
    expect_identical(predict(object, testx), ans)

    perm <- rev(seq_along(case$time))
    object$x <- object$x[perm, , drop = FALSE]
    object$y <- object$y[perm]
    object$censor <- object$censor[perm]
    expect_equal(predict(object, testx), ans, tolerance = 1e-12)
  }
})

test_that("survival reverse sums retain small late risk masses", {
  x <- matrix(c(0, rep(1, 5)), ncol = 1)
  bw <- 1 / (sd(x[, 1]) * sqrt(2 * 20 * log(10)))
  object <- revision_fit_object(x, 1:6, matrix(1), bw, c(0, rep(1, 5)))
  ans <- predict(object, matrix(0))
  expect_equal(ans$surv[, 1], c(.8, .6, .4, .2, 0), tolerance = 1e-12)
  expect_identical(ans$timepoints, 2:6)

  x <- matrix(c(0, 0, 100, 100), ncol = 1)
  object <- revision_fit_object(x, 1:4, matrix(1), .001, rep(1, 4))
  expect_equal(predict(object, matrix(0))$surv[, 1], c(.5, 0, 0, 0))
  expect_warning(ans <- predict(object, matrix(c(0, 1e6), ncol = 1)), "underflowed to zero")
  expect_equal(ans$surv[, 1], c(.5, 0, 0, 0))
  expect_true(all(is.na(ans$surv[, 2])))
})

test_that("invalid survival masses are diagnosed instead of replaced", {
  local_mocked_bindings(
    prediction_kernel = function(...) matrix(c(1, -2), 1, 2),
    .package = "orthoDr"
  )
  object <- revision_fit_object(matrix(1:2), 1:2, matrix(1), .5, c(1, 1))
  expect_error(predict(object, matrix(0)), "Invalid weighted survival")
})

test_that("public fitted models predict with unchanged fitted fields", {
  set.seed(7103)
  n <- 50L
  x <- sweep(matrix(rnorm(n * 3), n, 3), 2, c(1, 5, .3), "*")
  x <- sweep(x, 2, c(3, -9, 6), "+")
  y <- x[, 1] + rnorm(n)
  B <- matrix(c(1, 0, 0), 3, 1)
  fit <- orthoDr_reg(x, y, method = "sir", ndr = 1, B.initial = B,
                    bw = .6, keep.data = TRUE, maxitr = 5, ncore = 1)
  original <- unserialize(serialize(fit, NULL))
  tx <- x[1:2, , drop = FALSE]
  k <- revision_kernel_reference(x, tx, fit$B, fit$bw)
  expect_equal(predict(fit, tx)$pred, drop(k %*% y / rowSums(k)), tolerance = 1e-12)
  expect_identical(fit, original)

  time <- exp(.2 * x[, 1] + rnorm(n))
  event <- rep(c(1, 1, 0, 1, 0), length.out = n)
  forward <- orthoDr_surv(x, time, event, method = "forward", ndr = 1,
                         B.initial = B, bw = .6, keep.data = TRUE,
                         maxitr = 1, ncore = 1)
  expect_identical(forward$bw, .6)
  ans <- predict(forward, tx)
  expect_equal(dim(ans$surv), c(sum(event), nrow(tx)))
  expect_true(all(is.finite(ans$surv)))
})

test_that("all initializer candidates have exactly the requested dimension", {
  p <- 4L
  x <- matrix(seq_len(20 * p), 20, p)
  y <- seq_len(20)
  seen <- list()
  local_mocked_bindings(
    dr = function(formula, method) {
      order <- if (method == "save") c(2L, 1L, 3L, 4L) else seq_len(p)
      list(evectors = diag(p)[, order, drop = FALSE])
    },
    reg_init = function(method, Bc, ...) {
      seen[[length(seen) + 1L]] <<- Bc
      # The old two-column candidate would win when ndr exceeds two.
      if (ncol(Bc) == 2L) -10 else 1
    },
    .package = "orthoDr"
  )
  for (d in seq_len(p)) {
    seen <- list()
    B <- orthoDr:::initB(x, y, d, .4, "sir", 1L)
    expect_equal(dim(B), c(p, d))
    expect_equal(length(seen), if (d == 2L) 4L else 3L)
    expect_true(all(vapply(seen, function(candidate) identical(dim(candidate), c(p, d)), logical(1))))
  }
})

test_that("a dependent mixed initializer leaves the three valid starters", {
  calls <- 0L
  local_mocked_bindings(
    dr = function(...) list(evectors = diag(3)),
    reg_init = function(method, Bc, ...) { calls <<- calls + 1L; 1 },
    .package = "orthoDr"
  )
  B <- orthoDr:::initB(matrix(seq_len(60), 20, 3), 1:20, 2L, .4, "sir", 1L)
  expect_equal(calls, 3L)
  expect_equal(dim(B), c(3L, 2L))
})

test_that("real dr and pracma default initializers preserve dimensions", {
  set.seed(7104)
  x <- matrix(rnorm(100 * 4), 100, 4)
  y <- x[, 1] + x[, 2]^2 + rnorm(100)
  for (d in c(2L, 3L, 4L)) {
    B <- orthoDr:::initB(x, as.matrix(scale(y)), d, .5, "sir", 1L)
    expect_equal(dim(B), c(4L, d))
    expect_equal(crossprod(B), diag(d), tolerance = 1e-10)
  }
  fit <- orthoDr_reg(x, y, ndr = 3, method = "sir", maxitr = 5, ncore = 1)
  expect_equal(dim(fit$B), c(4L, 3L))
  expect_true(all(is.finite(fit$B)))
})

test_that("pSAVE bounds the central interval and preserves threshold multiplicity", {
  groups <- list()
  moments <- list()
  local_mocked_bindings(
    dr = function(formula, data, ...) {
      groups[[length(groups) + 1L]] <<- data$a
      M <- matrix(c(sum(data$a), 1, 1, sum(!data$a)), 2, 2)
      moments[[length(moments) + 1L]] <<- M
      list(M = M)
    },
    .package = "orthoDr"
  )
  for (n in c(70L, 99L, 100L, 101L, 200L, 201L)) {
    x <- cbind(seq_len(n), sin(seq_len(n)))
    a <- rep(seq_len(ceiling(n / 2)), each = 2)[seq_len(n)]
    groups <- moments <- list()
    B <- pSAVE(x, a, seq_len(n), ndr = 1)
    indices <- seq.int(max(1L, floor(n / 2 - 50)), min(n, floor(n / 2 + 50)))
    thresholds <- sort(a)[indices]
    if (n != 70L)
      expect_identical(thresholds, sort(a)[(n / 2 - 50):(n / 2 + 50)])
    admissible <- vapply(thresholds, function(cut) {
      group <- as.numeric(a <= cut)
      all(vapply(unique(group), function(g) {
        xx <- x[group == g, , drop = FALSE]
        nrow(xx) >= 3L && qr(sweep(xx, 2L, colMeans(xx), "-"))$rank == 2L
      }, logical(1)))
    }, logical(1))
    thresholds <- thresholds[admissible]
    expect_identical(groups, lapply(thresholds, function(cut) as.numeric(a <= cut)))
    expected <- eigen(Reduce("+", moments), symmetric = FALSE)$vectors[, 1, drop = FALSE]
    expect_equal(tcrossprod(B), tcrossprod(expected), tolerance = 1e-12)
  }
})

test_that("pSAVE retains real grouped SAVE behavior on feasible inputs", {
  set.seed(7105)
  n <- 200L
  x <- matrix(rnorm(n * 3), n, 3)
  y <- x[, 1] + rnorm(n)
  a <- runif(n)
  for (dose in list(a, round(a, 1), rep(1, n))) {
    B <- pSAVE(x, dose, y, ndr = 1)
    expect_equal(dim(B), c(3L, 1L))
    expect_true(all(is.finite(B)))
  }
  B <- pSAVE(x[1:70, , drop = FALSE], rep(1, 70), y[1:70], ndr = 1)
  expect_equal(dim(B), c(3L, 1L))
  expect_true(all(is.finite(B)))
})

test_that("hMave accepts documented vectors and preserves matrix results", {
  set.seed(7106)
  n <- 24L
  x <- matrix(rnorm(n * 2), n, 2)
  y <- rexp(n, exp(.2 * x[, 1]))
  event <- rep(c(TRUE, TRUE, FALSE), length.out = n)
  matrices <- hMave(x, as.matrix(y), as.matrix(as.numeric(event)), 1)
  expect_equal(hMave(x, y, event, 1), matrices, tolerance = 1e-12)
  expect_equal(hMave(x, y, as.numeric(event), 1), matrices, tolerance = 1e-12)
  expect_error(hMave(x, y[-1], event, 1), "one value per observation")
  expect_error(hMave(x, y, event[-1], 1), "one value per observation")
  expect_error(hMave(x, cbind(y, y), event, 1), "one-column matrix")
  expect_error(hMave(x, y, cbind(event, event), 1), "one-column matrix")
})

test_that("distance vector coercion and abbreviated dispatch preserve metrics", {
  s1 <- c(1, 0, 0)
  s2 <- c(sqrt(3) / 2, .5, 0)
  x <- rbind(diag(3), -diag(3))
  prefixes <- c(dist = "d", trace = "tr", canonical = "can", sine = "sin")
  for (type in names(prefixes)) {
    reference <- distance(as.matrix(s1), as.matrix(s2), type = type, x = x)
    expect_equal(distance(s1, s2, type = type, x = x), reference)
    expect_equal(distance(s1, s2, type = prefixes[[type]], x = x), reference)
  }
  expect_equal(distance(s1, s2, "dist"), sqrt(.5), tolerance = 1e-12)
  expect_equal(distance(s1, s2, "trace"), .75, tolerance = 1e-12)
  expect_equal(distance(s1, s2, "canonical", x), sqrt(.75), tolerance = 1e-12)
  expect_equal(distance(s1, s2, "sine"), .25, tolerance = 1e-12)
  expect_warning(a <- distance(matrix(c(1, 0)), diag(2), "sine"), "Dimension d")
  expect_warning(b <- distance(diag(2), matrix(c(1, 0)), "sine"), "Dimension d")
  expect_equal(a, 0)
  expect_equal(b, 1)
})
