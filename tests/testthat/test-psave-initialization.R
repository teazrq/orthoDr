# Independent Cholesky local coordinates and a raw coefficient eigenproblem.
# Only the separately tested slice partition is shared with the implementation.
psave_coordinate_reference <- function(x, a, r, ndr = 1, nslices0 = 2,
                                        transport_fallback = FALSE) {
  n <- nrow(x)
  p <- ncol(x)
  cuts <- sort(a)[seq.int(max(1L, floor(n / 2 - 50)),
                          min(n, floor(n / 2 + 50)))]
  reference <- raw <- matrix(0, p, p)
  fallback <- 0
  retained <- numeric()
  for (cut in cuts) {
    group <- a <= cut
    indices <- lapply(unique(group), function(g) which(group == g))
    if (!all(vapply(indices, function(i) {
      length(i) >= p + 1L && qr(sweep(x[i, , drop = FALSE], 2L,
                                      colMeans(x[i, , drop = FALSE]), "-"))$rank == p
    }, logical(1)))) next
    retained <- c(retained, cut)
    for (i in indices) {
      xx <- x[i, , drop = FALSE]
      centered <- sweep(xx, 2L, colMeans(xx), "-")
      covariance <- crossprod(centered) / length(i)
      inverse <- solve(t(chol(covariance)))
      slices <- orthoDr:::.psave_slices(r[i], nslices0)
      reference <- reference + length(i) / n * covariance
      for (j in seq_len(slices$nslices)) {
        size <- slices$slice.sizes[j]
        if (size < 3L) {
          fallback <- fallback + size / n
          if (transport_fallback)
            raw <- raw + size / n * crossprod(inverse)
          next
        }
        subset <- xx[slices$slice.indicator == j, , drop = FALSE]
        subset <- sweep(subset, 2L, colMeans(subset), "-")
        within <- crossprod(subset) / size
        contrast <- diag(p) - inverse %*% within %*% t(inverse)
        raw <- raw + size / n * t(inverse) %*% contrast %*% contrast %*% inverse
      }
    }
  }
  reference <- reference / length(retained)
  raw <- raw / length(retained)
  solution <- eigen(raw %*% reference, symmetric = FALSE)
  order <- order(Re(solution$values), decreasing = TRUE)
  directions <- Re(solution$vectors[, order[seq_len(ndr)], drop = FALSE])
  directions <- sweep(directions, 2L, sqrt(colSums(directions^2)), "/")
  for (j in seq_len(ndr))
    if (directions[which.max(abs(directions[, j])), j] < 0)
      directions[, j] <- -directions[, j]
  list(B = directions, values = Re(solution$values[order]), covariance = reference,
       raw.moment = raw, fallback.mass = fallback / length(retained), thresholds = retained)
}

psave_projector <- function(b) tcrossprod(qr.Q(qr(b)))

test_that("pSAVE validates real finite inputs and dimensions before initialization", {
  set.seed(261030)
  x <- matrix(rnorm(60), 20, 3)
  a <- seq_len(20)
  r <- x[, 1]
  expect_error(pSAVE(as.data.frame(x), a, r), "X must be a matrix")
  expect_error(pSAVE(x + 1i, a, r), "real")
  expect_error(pSAVE(x[-1, ], a, r), "Number of observations do not match")
  bad <- x
  bad[1, 1] <- Inf
  expect_error(pSAVE(bad, a, r), "finite numeric values")
  expect_error(pSAVE(matrix(numeric(), 20, 0), a, r), "at least two rows and one column")
  expect_error(pSAVE(matrix(1, 1, 1), 1, 1, ndr = 1), "at least two rows")
  for (value in list(as.character(a), factor(a), a + 1i, replace(a, 1, NA),
                     replace(a, 1, Inf), matrix(a, ncol = 1)))
    expect_error(pSAVE(x, value, r), "a must be a finite numeric vector")
  for (value in list(as.character(r), r + 1i, replace(r, 1, NaN),
                     matrix(r, ncol = 1)))
    expect_error(pSAVE(x, a, value), "r must be a finite numeric vector")
  for (value in list(0, 4, 1.5, NA_real_, Inf, TRUE, c(1, 2), matrix(1)))
    expect_error(pSAVE(x, a, r, ndr = value), "ndr must be an integer")
  for (value in list(0, -1, 1.5, NA_real_, Inf, TRUE, c(1, 2), matrix(2)))
    expect_error(pSAVE(x, a, r, nslices0 = value), "nslices0 must be a positive integer")
})

test_that("the slicing adapter preserves successful dependency results", {
  set.seed(261026)
  old_success <- old_failure <- preserved <- corrected <- 0L
  valid_contract <- TRUE
  for (n in 2:30) for (h in 1:6) {
    rewards <- list(seq_len(n), rep(0:2, length.out = n),
                    sample(0:3, n, replace = TRUE), c(rep(0, n - 1), 1),
                    rnorm(n))
    for (y in rewards) {
      old <- tryCatch(dr::dr.slices(y, h), error = identity)
      new <- orthoDr:::.psave_slices(y, h)
      valid_contract <- valid_contract &&
        length(new$slice.indicator) == n &&
        length(new$slice.sizes) == new$nslices &&
        sum(new$slice.sizes) == n &&
        all(tabulate(new$slice.indicator, nbins = new$nslices) == new$slice.sizes)
      if (inherits(old, "error")) {
        old_failure <- old_failure + 1L
        corrected <- corrected + as.integer(new$nslices == 1L)
      } else {
        old_success <- old_success + 1L
        preserved <- preserved + as.integer(identical(new, old))
      }
    }
  }
  expect_true(valid_contract)
  expect_gt(old_success, 0L)
  expect_gt(old_failure, 0L)
  expect_identical(preserved, old_success)
  expect_identical(corrected, old_failure)
  repaired <- orthoDr:::.psave_slices(c(rep(0, 18), 1, 2), 2)
  expect_identical(repaired, list(slice.indicator = rep(1, 20), nslices = 1L,
                                  slice.sizes = 20L))
})

test_that("corrected pSAVE matches independent raw-coordinate moments with dose ties", {
  set.seed(261028)
  for (n in c(150L, 200L)) {
    x <- matrix(rnorm(n * 3), n, 3)
    r <- x[, 1] + rnorm(n)
    a <- runif(n)
    for (dose in list(a, round(a, 1), rep(1, n))) {
      actual <- pSAVE(x, dose, r, ndr = 1)
      expected <- psave_coordinate_reference(x, dose, r)
      expect_gt(expected$values[1] - expected$values[2], 1e-5)
      expect_equal(actual, expected$B, tolerance = 2e-11)
      expect_equal(colSums(actual^2), 1, tolerance = 1e-14)
      expect_gt(actual[which.max(abs(actual[, 1])), 1], 0)
    }
  }
})

test_that("singleton and uninformative groups preserve weights in common coordinates", {
  set.seed(261031)
  x <- matrix(rnorm(60), 20, 3)
  a <- rep(1, 20)
  r <- c(rep(0, 19), 1)
  expect_identical(orthoDr:::.psave_slices(r, 2)$slice.sizes, c(19, 1))
  expect_equal(pSAVE(x, a, r, ndr = 1), psave_coordinate_reference(x, a, r)$B,
               tolerance = 2e-11)

  # The first ten thresholds have only constant-outcome groups. The final ten
  # have a single group with outcome contrasts. Preserve all twenty weights.
  a <- rep(0:1, each = 10)
  r <- a
  expect_equal(pSAVE(x, a, r, ndr = 1), psave_coordinate_reference(x, a, r)$B,
               tolerance = 2e-11)

  # A singleton-only group contributes identity alongside an informative group.
  x <- matrix(rnorm(120), 40, 3)
  a <- rep(0:1, each = 20)
  r <- c(seq_len(20), rep(0:1, each = 10))
  expect_equal(pSAVE(x, a, r, ndr = 1, nslices0 = 20),
                psave_coordinate_reference(x, a, r, nslices0 = 20)$B,
                tolerance = 2e-11)
})

test_that("local SAVE matches dependency slice covariances and separates tiny fallback", {
  set.seed(261033)
  x <- matrix(rnorm(60), 20, 3) %*% diag(c(1, 4, .7))
  r <- c(rep(0, 19), 1)
  local <- orthoDr:::.psave_group_moment(x, r, 2)
  expected <- dr::dr(r ~ x, data = list(x = x, a = rep(1, 20), r = r),
                     group = ~a, nslices = 2, numdir = 3, method = "save",
                     slice.function = orthoDr:::.psave_slices)
  expect_equal(local$moment + local$fallback.mass * diag(3), expected$M,
               tolerance = 2e-12)
  expect_equal(local$fallback.mass, 1 / 20)
  expect_equal(local$covariance, crossprod(sweep(x, 2L, colMeans(x), "-")) / 20,
               tolerance = 1e-14)
  tiny <- orthoDr:::.psave_group_moment(x, seq_len(20), 20)
  expect_identical(tiny$moment, matrix(0, 3, 3))
  expect_equal(tiny$fallback.mass, 1)
})

test_that("equal-covariance correction recovers an anisotropic mixed raw direction", {
  covariance <- matrix(c(2, 1.2, .3, 1.2, 16, 1, .3, 1, .7), 3, 3)
  eig <- eigen(covariance, symmetric = TRUE)
  root <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  beta <- c(1, 1, 0) / sqrt(2)
  q <- qr.Q(qr(cbind(root %*% beta, c(0, 0, 1), c(1, 0, 0))))
  blocks <- lapply(c(.25, 1.75), function(v) {
    z <- diag(sqrt(3 * c(v, 1, 1)))
    rbind(z, -z) %*% t(q) %*% root
  })
  x <- do.call(rbind, blocks)
  r <- rep(0:1, each = 6)
  actual <- pSAVE(x, rep(1, 12), r, ndr = 1)
  expect_equal(actual[, 1], beta, tolerance = 5e-14)
  local <- orthoDr:::.psave_group_moment(x, r, 2)
  reduced <- local$inverse.sqrt %*% eigen(local$moment, symmetric = TRUE)$vectors[, 1]
  expect_equal(psave_projector(actual), psave_projector(matrix(reduced, 3, 1)),
               tolerance = 5e-14)
})

test_that("anisotropic correlated Gaussian direction is recovered approximately", {
  set.seed(716)
  scales <- c(1, 5, 2, .5)
  covariance <- outer(scales, scales) * .6^abs(outer(seq_len(4), seq_len(4), "-"))
  x <- matrix(rnorm(4000 * 4), 4000, 4) %*% chol(covariance)
  beta <- c(1, 1, 0, 0) / sqrt(2)
  r <- as.numeric(x %*% beta) ^ 2 / as.numeric(t(beta) %*% covariance %*% beta) +
    rnorm(nrow(x), sd = .02)
  actual <- pSAVE(x, rep(1, nrow(x)), r, ndr = 1)
  expect_lt(sqrt(sum((psave_projector(actual) - tcrossprod(beta))^2)), .12)
})

test_that("a complete tied eigenspace is compared through its projector", {
  covariance <- matrix(c(2, 1.2, .3, 1.2, 16, 1, .3, 1, .7), 3, 3)
  eig <- eigen(covariance, symmetric = TRUE)
  root <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  blocks <- lapply(c(.25, 1.75), function(v) {
    z <- diag(sqrt(3 * c(v, v, 1)))
    rbind(z, -z) %*% root
  })
  x <- do.call(rbind, blocks)
  r <- rep(0:1, each = 6)
  expected <- solve(root, diag(3)[, 1:2])
  actual <- pSAVE(x, rep(1, 12), r, ndr = 2)
  expect_equal(psave_projector(actual), psave_projector(expected), tolerance = 5e-14)
  expect_equal(psave_projector(pSAVE(x[12:1, ], rep(1, 12), r[12:1], ndr = 2)),
               psave_projector(expected), tolerance = 5e-14)
})

test_that("unequal group moments transform affinely in raw selected coordinates", {
  set.seed(119)
  n <- 181
  a <- round(runif(n, 0, 4), 1)
  x <- matrix(rnorm(n * 4), n, 4)
  transform <- diag(c(1, 3, .7, 2))
  transform[1, 2] <- .8
  high <- a > 2
  x[high, ] <- sweep(x[high, ] %*% transform, 2L, 0:3, "+")
  r <- (x[, 1] + .4 * x[, 2])^2 + .3 * x[, 3] + rnorm(n)
  expected <- psave_coordinate_reference(x, a, r, ndr = 2)
  expect_gt(expected$values[2] - expected$values[3], 1e-4)
  actual <- pSAVE(x, a, r, ndr = 2)
  expect_equal(actual, expected$B, tolerance = 3e-11)
  change <- matrix(c(2, 0, .2, 0, .3, .7, 0, -.3,
                      0, .4, 1.4, 0, -.2, 0, .1, 1.1), 4, 4)
  shifted <- sweep(x %*% change, 2L, c(5, -2, 1, 8), "+")
  transformed <- pSAVE(shifted, a, r, ndr = 2)
  expect_equal(psave_projector(transformed), psave_projector(solve(change, actual)),
               tolerance = 3e-11)
})

test_that("mixed tiny slices are direction neutral in common coordinates", {
  set.seed(88)
  x <- matrix(rnorm(300 * 3), 300, 3)
  x[1:250, ] <- x[1:250, ] %*% matrix(c(1, 0, 0, 3, 6, 0, 0, 0, .4), 3, 3)
  x[251:300, ] <- sweep(x[251:300, ], 2L, c(2, -1, 4), "+")
  a <- rep(0:1, c(250, 50))
  r <- c(1:250, rep(0:1, c(20, 30)))
  expected <- psave_coordinate_reference(x, a, r, nslices0 = 300)
  actual <- pSAVE(x, a, r, ndr = 1, nslices0 = 300)
  expect_equal(expected$fallback.mass, 250 / 300)
  expect_equal(psave_projector(actual), psave_projector(expected$B), tolerance = 3e-11)
  polluted <- psave_coordinate_reference(x, a, r, nslices0 = 300,
                                          transport_fallback = TRUE)
  expect_gt(sqrt(sum((psave_projector(actual) - psave_projector(polluted$B))^2)), .2)
  original <- orthoDr:::.psave_grouped_moment
  local_mocked_bindings(
    .psave_grouped_moment = function(...) {
      local <- original(...)
      local$fallback.mass <- local$fallback.mass + 100
      local
    }, .package = "orthoDr"
  )
  expect_identical(pSAVE(x, a, r, ndr = 1, nslices0 = 300), actual)
})

test_that("central thresholds retain repeated weights and odd-size truncation", {
  set.seed(133)
  x <- matrix(rnorm(151 * 3), 151, 3)
  a <- round(runif(151), 1)
  r <- x[, 1] + rnorm(151)
  observed <- list()
  original <- orthoDr:::.psave_grouped_moment
  local_mocked_bindings(
    .psave_grouped_moment = function(x, y, group, nslices) {
      observed[[length(observed) + 1L]] <<- group
      original(x, y, group, nslices)
    }, .package = "orthoDr"
  )
  pSAVE(x, a, r, ndr = 1)
  cuts <- sort(a)[25:125]
  expect_length(observed, 101L)
  expect_lt(length(unique(cuts)), 101L)
  expect_identical(observed, lapply(cuts, function(cut) as.numeric(a <= cut)))
})

test_that("pSAVE rejects structural zero and identity slice information", {
  set.seed(261027)
  expect_error(pSAVE(rbind(c(0, 0), c(1, 0), c(0, 1)), 1:3, c(0, 0, 1),
                     ndr = 1), "no directional slice information.*B.initial")
  x <- matrix(rnorm(60), 20, 3)
  expect_error(pSAVE(x, rep(1, 20), 1:20, ndr = 1, nslices0 = 20),
               "no directional slice information.*B.initial")
  expect_error(pSAVE(x, rep(1, 20), rep(0, 20), ndr = 1),
               "no outcome-slice contrasts.*B.initial")
  x <- matrix(rnorm(900), 300, 3)
  expect_error(pSAVE(x, c(rep(0, 250), rep(1, 50)),
                     c(rep(0, 250), 1:50), ndr = 1, nslices0 = 300),
               "no directional slice information.*B.initial")
})

test_that("pSAVE distinguishes global rank failure from no admissible thresholds", {
  set.seed(261032)
  x <- matrix(rnorm(120), 40, 3)
  x[, 3] <- x[, 1] + x[, 2]
  expect_error(pSAVE(x, 1:40, rnorm(40), ndr = 1),
               "full centred predictor rank.*B.initial")
  expect_error(pSAVE(matrix(rnorm(30), 5, 6), 1:5, rnorm(5), ndr = 1),
               "full centred predictor rank.*B.initial")
  # Full global rank, but every central cutoff splits 104 rows into groups
  # that cannot both contain the required 55 rows. The maximum is not central.
  x <- matrix(rnorm(104 * 54), 104, 54)
  expect_error(pSAVE(x, 1:104, rnorm(104), ndr = 1),
               "no admissible dose threshold.*B.initial")
})

test_that("each actual dose group must have full centred predictor rank", {
  x <- cbind(1:12, c(2 * (1:6), sin(7:12)))
  a <- 1:12
  observed <- list()
  original <- orthoDr:::.psave_grouped_moment
  local_mocked_bindings(
    .psave_grouped_moment = function(x, y, group, nslices) {
      observed[[length(observed) + 1L]] <<- group
      original(x, y, group, nslices)
    },
    .package = "orthoDr"
  )
  B <- pSAVE(x, a, 1:12, ndr = 1)
  expect_true(all(is.finite(B)))
  expect_identical(observed, lapply(c(7, 8, 9, 12), function(cut) {
    as.numeric(a <= cut)
  }))
})

test_that("one predictor returns its exact subspace only after validation", {
  x <- matrix(1:3, 3, 1)
  expect_identical(pSAVE(x, 1:3, rep(1, 3), ndr = 1), matrix(1, 1, 1))
  expect_error(pSAVE(x, 1:3, 1:3, ndr = 2), "ndr must be an integer")
  expect_error(pSAVE(x, 1:3, 1:3, ndr = 1, nslices0 = 0), "positive integer")
  expect_error(pSAVE(x, 1:3, c(1, NA, 3), ndr = 1), "finite numeric vector")
  expect_error(pSAVE(matrix(1, 3, 1), 1:3, 1:3, ndr = 1),
               "full centred predictor rank")
})

test_that("small-sample automatic direct and pseudo-direct fits use repaired pSAVE", {
  set.seed(127)
  x <- matrix(runif(60 * 4, -1, 1), 60, 4)
  a <- runif(60, 0, 2)
  target <- 1 + .4 * x[, 1] - .3 * x[, 2]
  r <- 3 + x[, 1] - 2 * (a - target)^2 + rnorm(60, sd = .2)
  for (n in c(45L, 60L)) {
    xx <- x[seq_len(n), , drop = FALSE]
    aa <- a[seq_len(n)]
    rr <- r[seq_len(n)]
    initial <- pSAVE(xx, aa, rr, ndr = 2)
    for (method in c("direct", "pseudo_direct")) {
      args <- list(x = xx, a = aa, r = rr, ndr = 2, method = method,
                   K = 5, keep.data = TRUE, maxitr = 1, ncore = 1)
      automatic <- do.call(orthoDr_pdose, args)
      supplied <- do.call(orthoDr_pdose, c(args, list(B.initial = initial)))
      expect_true(all(is.finite(automatic$B)))
      expect_identical(automatic$B, supplied$B)
      expect_identical(predict(automatic, xx)$pred, predict(supplied, xx)$pred)
    }
  }
})

test_that("supplied starts bypass full predictor rank requirements in both methods", {
  set.seed(261029)
  x <- matrix(rnorm(24), 6, 4)
  x[, 4] <- x[, 1] + x[, 2]
  wide <- matrix(rnorm(30), 5, 6)
  for (xx in list(x, wide)) {
    a <- runif(nrow(xx))
    r <- rnorm(nrow(xx))
    initial <- matrix(c(1, rep(0, ncol(xx) - 1L)), ncol(xx), 1)
    expect_error(pSAVE(xx, a, r, ndr = 1), "full centred predictor rank")
    for (method in c("direct", "pseudo_direct")) {
      fit <- orthoDr_pdose(xx, a, r, ndr = 1, B.initial = initial,
                          method = method, K = 3, maxitr = 1, ncore = 1)
      expect_true(all(is.finite(fit$B)))
    }
  }
})

test_that("one-predictor automatic fitting matches the exact supplied basis", {
  x <- matrix(seq(-1, 1, length.out = 12), 12, 1)
  a <- seq(.1, 1.9, length.out = 12)
  r <- x[, 1] + a
  for (method in c("direct", "pseudo_direct")) {
    args <- list(x = x, a = a, r = r, ndr = 1, method = method,
                 K = 3, keep.data = TRUE, maxitr = 1, ncore = 1)
    automatic <- do.call(orthoDr_pdose, args)
    supplied <- do.call(orthoDr_pdose, c(args, list(B.initial = matrix(1, 1, 1))))
    expect_identical(automatic$B, supplied$B)
    expect_identical(predict(automatic, x)$pred, predict(supplied, x)$pred)
  }
})
