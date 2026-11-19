test_that("direct-dose fitting and its basis gradient use predictor differences", {
  i <- 0:31
  x <- cbind((11 * i + 4) %% 33 - 16, (7 * i^2 + 13 * i + 2) %% 33 - 16)
  query <- cbind((17 * i + 3) %% 33 - 16, (19 * i + 7) %% 33 - 16)
  a <- ((23 * i + 5) %% 37) / 36
  r <- -(a - .5 - .3 * sin(rowSums(x) / 8))^2
  fit <- function(z) orthoDr_pdose(z, a, r, ndr = 1, B.initial = matrix(c(1, 0)),
                                  intercept = TRUE, keep.data = TRUE, ncore = 1)
  base <- fit(x)
  evaluate <- function(z) orthoDr:::pdose_direct_evaluate(
    base$B, z, a, r, base$bw, base$W, .1, intercept = TRUE, alpha = base$intercept)
  state <- evaluate(x)
  expect_true(state$valid)
  for (offset in c(1e9, 2^52)) {
    shifted <- fit(x + offset)
    moved_state <- evaluate(x + offset)
    expect_true(moved_state$valid)
    expect_equal(moved_state$gB, state$gB, tolerance = 1e-12)
    expect_equal(shifted$B, base$B, tolerance = 1e-12)
    expect_equal(shifted$W, base$W, tolerance = 1e-12)
    expect_identical(shifted$stop.reason, base$stop.reason)
    expect_equal(predict(shifted, query + offset)$pred, predict(base, query)$pred,
                 tolerance = 1e-12)
  }
})

test_that("dose prediction resolves small index differences after an exact shift", {
  grid <- expand.grid(u = c(-8, 0, 8), z = -1:1, a = seq(0, 1, length.out = 9))
  x <- cbind(grid$u, grid$z - grid$u)
  a <- grid$a
  reward <- -(a - (.5 + .4 * grid$z))^2
  query <- cbind(0, -1:1)
  for (method in c("direct", "pseudo_direct")) for (intercept in c(FALSE, TRUE)) {
    fit <- function(z) orthoDr_pdose(z, a, reward, ndr = 1, bw = .6, K = 9,
      method = method, B.initial = matrix(c(1, 1)), intercept = intercept,
      keep.data = TRUE, ncore = 1)
    base <- fit(x)
    shifted <- fit(sweep(x, 2L, c(2^55, 0), "+"))
    train_index <- drop(x %*% shifted$B)
    query_index <- drop(query %*% shifted$B)
    K <- exp(-.5 * (outer(query_index, train_index, "-") / (sd(train_index) * shifted$bw))^2)
    expected <- drop(K %*% shifted$W)
    if (intercept) expected <- expected + shifted$intercept
    if (intercept && method == "direct")
      expected <- expected - sum(shifted$kernel.mean * shifted$W)
    actual <- drop(predict(shifted, sweep(query, 2L, c(2^55, 0), "+"))$pred)
    expect_gt(diff(range(expected)), .1)
    expect_equal(actual, expected, tolerance = 1e-12)
    expect_equal(actual, drop(predict(base, query)$pred), tolerance = 1e-12)
  }
})

test_that("a large observation with a zero contrast preserves the dose kernel", {
  s <- rep(-1:1, 12)
  x <- rbind(c(2^55, 2^55), cbind(s, 0))
  b <- matrix(c(1, -1) / sqrt(2))
  query <- cbind(-1:1, 0)
  w <- c(0, -s / 12)
  z <- c(0, s) / sqrt(2)
  expected <- exp(-.5 * (outer((-1:1) / sqrt(2), z, "-") / (.5 * sd(z)))^2) %*% w
  for (order in list(seq_len(nrow(x)), c(2:nrow(x), 1L))) {
    actual <- orthoDr:::dosepred(b, x[order, ], query, .5, w[order])
    expect_equal(drop(actual), drop(expected), tolerance = 1e-12)
  }
})
