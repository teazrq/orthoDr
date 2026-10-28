psave_legacy_reference <- function(x, a, r, ndr = 1, nslices0 = 2) {
  n <- nrow(x)
  p <- ncol(x)
  cuts <- sort(a)[seq.int(max(1L, floor(n / 2 - 50)),
                          min(n, floor(n / 2 + 50)))]
  total <- matrix(0, p, p)
  for (cut in cuts) {
    fit <- dr::dr(r ~ x, data = list(x = x, a = as.numeric(a <= cut), r = r),
                  group = ~a, nslices = nslices0, chi2approx = "bx",
                  numdir = p, method = "save")
    total <- total + fit$M
  }
  eigen(total, symmetric = FALSE)$vectors[, seq_len(ndr), drop = FALSE]
}

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

test_that("feasible pSAVE retains historical tied and one-group moments exactly", {
  set.seed(261028)
  for (n in c(150L, 200L)) {
    x <- matrix(rnorm(n * 3), n, 3)
    r <- x[, 1] + rnorm(n)
    a <- runif(n)
    for (dose in list(a, round(a, 1), rep(1, n))) {
      actual <- pSAVE(x, dose, r, ndr = 1)
      expected <- psave_legacy_reference(x, dose, r)
      expect_identical(actual, expected)
    }
  }
})

test_that("singleton and uninformative groups retain their historical contributions", {
  set.seed(261031)
  x <- matrix(rnorm(60), 20, 3)
  a <- rep(1, 20)
  r <- c(rep(0, 19), 1)
  expect_identical(orthoDr:::.psave_slices(r, 2)$slice.sizes, c(19, 1))
  expect_identical(pSAVE(x, a, r, ndr = 1), psave_legacy_reference(x, a, r))

  # The first ten thresholds have only constant-outcome groups. The final ten
  # have a single group with outcome contrasts. Preserve all twenty weights.
  a <- rep(0:1, each = 10)
  r <- a
  expect_identical(pSAVE(x, a, r, ndr = 1), psave_legacy_reference(x, a, r))

  # A singleton-only group contributes identity alongside an informative group.
  x <- matrix(rnorm(120), 40, 3)
  a <- rep(0:1, each = 20)
  r <- c(seq_len(20), rep(0:1, each = 10))
  expect_identical(pSAVE(x, a, r, ndr = 1, nslices0 = 20),
                   psave_legacy_reference(x, a, r, nslices0 = 20))
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
  local_mocked_bindings(
    dr = function(formula, data, ...) {
      observed[[length(observed) + 1L]] <<- data$a
      list(M = matrix(c(sum(data$a), 1, 1, sum(!data$a)), 2, 2))
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
