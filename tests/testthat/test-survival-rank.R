test_that("survival fitting detects constant and dependent predictors in any order", {
  i <- 0:31
  a <- i %% 7 - 3
  b <- (5*i) %% 11 - 5
  x <- cbind(a, b, a+b)
  t <- (11*i) %% 37 + 1
  e <- i %% 3 != 1
  for (j in seq_len(3)) {
    k <- setdiff(1:3, j)
    for (perm in list(c(j,k), c(j,rev(k)))) {
      expect_error(CP_SIR(x[, perm], t, e), "rank")
      expect_error(orthoDr_surv(x[, perm], t, e, ndr=1), "rank")
    }
  }
  expect_error(orthoDr_surv(cbind(a,a), t, e, ndr=1,
                            B.initial=matrix(c(1,-1))), "projected")
  expect_error(orthoDr_surv(cbind(a, .1), t, e, ndr=1,
                            B.initial=matrix(c(1, 0))), "Constant predictor")
  expect_error(CP_SIR(cbind(a, .1), t, e), "Constant predictor")
  units <- c(1, 1e-12)
  orig <- CP_SIR(x[, 1:2], t, e)$vectors[, 1]
  scaled <- CP_SIR(sweep(x[, 1:2], 2, units, "*"), t, e)$vectors[, 1]
  expect_lt(distance(orig, scaled * units), 1e-11)
})

test_that("projector distances preserve equivalence under basis rescaling", {
  b <- diag(3)[, 1:2]
  other <- sweep(b, 2, c(1, 1e-25), "*")
  expect_equal(distance(b, other), 0)
  expect_equal(distance(b, other, "trace"), 1)
})

test_that("raw directions preserve the fitted line or diagnose incompatible units", {
  set.seed(46)
  x <- matrix(rnorm(128), 64, 2)
  t <- runif(64)
  e <- t < .7
  a <- runif(64)
  r <- runif(64)
  for (fit in list(function(z) CP_SIR(z, t, e)$vectors,
                   function(z) orthoDr_surv(z, t, e, ndr=2, maxitr=0,
                     B.initial=cbind(c(1,1), c(1,-1)))$B,
                   function(z) pSAVE(z, a, r, ndr=2))) {
    expect_error(fit(sweep(x, 2, c(1e-180, 1e170), "*")), "Rescale x")
    expect_lt(distance(fit(x)[,1], fit(x * 1e-312)[,1]), 1e-8)
  }
})
