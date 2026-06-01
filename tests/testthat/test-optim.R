# Tests for ortho_optim (general solver on Stiefel manifold)
# Strategy: eigenvalue problems with known analytical solutions

test_that("ortho_optim solves eigenvalue problem with analytical gradient", {
  set.seed(42)
  n <- 20
  k <- 3

  # Symmetric PD matrix: A = V D V'
  A <- matrix(rnorm(n * n), n, n)
  A <- crossprod(A)

  # Truth: top-k eigenvectors minimize -0.5 * tr(B'AB)
  eig <- eigen(A, symmetric = TRUE)
  B_true <- eig$vectors[, 1:k]
  f_true <- -0.5 * sum(eig$values[1:k])

  # Random orthogonal initial value
  B0 <- qr.Q(qr(matrix(rnorm(n * k), n, k)))

  # Objective and analytical gradient
  fn <- function(B) -0.5 * sum(diag(crossprod(B, A %*% B)))
  gr <- function(B) -A %*% B

  fit <- ortho_optim(B0, fn, gr, maxitr = 500, verbose = FALSE)

  # Objective should match truth
  expect_true(fit$converge)
  expect_equal(fit$fn, f_true, tolerance = 1e-4)

  # B should be orthonormal: diag(B'B) = 1, off-diag ~ 0
  BtB <- crossprod(fit$B)
  expect_equal(diag(as.matrix(BtB)), rep(1, k), tolerance = 1e-6)
  expect_equal(max(abs(BtB - diag(k))), 0, tolerance = 1e-8)

  # Subspace should match true eigenvectors via projection (n x n)
  proj_true <- tcrossprod(B_true)
  proj_fit  <- tcrossprod(fit$B)
  expect_equal(as.numeric(proj_fit), as.numeric(proj_true), tolerance = 1e-3)

  # Method indicator
  expect_equal(fit$method, "true gradient")
})

test_that("ortho_optim solves eigenvalue problem with numerical gradient", {
  set.seed(142)
  n <- 15
  k <- 2

  A <- matrix(rnorm(n * n), n, n)
  A <- crossprod(A)

  eig <- eigen(A, symmetric = TRUE)
  B_true <- eig$vectors[, 1:k]
  f_true <- -0.5 * sum(eig$values[1:k])

  B0 <- qr.Q(qr(matrix(rnorm(n * k), n, k)))

  fn <- function(B) -0.5 * sum(diag(crossprod(B, A %*% B)))

  fit <- ortho_optim(B0, fn, grad = NULL, maxitr = 500, verbose = FALSE)

  expect_true(fit$converge)
  expect_equal(fit$fn, f_true, tolerance = 1e-3)
  expect_equal(fit$method, "approx. gradient")

  # Subspace recovery
  proj_true <- tcrossprod(B_true)
  proj_fit  <- tcrossprod(fit$B)
  expect_equal(as.numeric(proj_fit), as.numeric(proj_true), tolerance = 1e-2)
})

test_that("ortho_optim maximize flag works", {
  set.seed(242)
  n <- 10
  k <- 2

  A <- matrix(rnorm(n * n), n, n)
  A <- crossprod(A)

  eig <- eigen(A, symmetric = TRUE)
  f_true <- 0.5 * sum(eig$values[1:k]) # maximized value (positive)

  B0 <- qr.Q(qr(matrix(rnorm(n * k), n, k)))

  # Maximize 0.5 * tr(B'AB)
  fn <- function(B) 0.5 * sum(diag(crossprod(B, A %*% B)))
  gr <- function(B) A %*% B

  fit <- ortho_optim(B0, fn, gr, maximize = TRUE, maxitr = 500, verbose = FALSE)

  expect_true(fit$converge)
  # fn should be the maximized value (not negated)
  expect_equal(fit$fn, f_true, tolerance = 1e-4)
})

test_that("ortho_optim works for ndr = 1 (single direction)", {
  set.seed(342)
  n <- 12

  A <- matrix(rnorm(n * n), n, n)
  A <- crossprod(A)

  eig <- eigen(A, symmetric = TRUE)
  v_true <- eig$vectors[, 1]
  f_true <- -0.5 * eig$values[1]

  B0 <- matrix(rnorm(n), n, 1)
  B0 <- B0 / sqrt(drop(crossprod(B0)))

  fn <- function(B) -0.5 * sum(diag(crossprod(B, A %*% B)))
  gr <- function(B) -A %*% B

  fit <- ortho_optim(B0, fn, gr, maxitr = 500, verbose = FALSE)

  expect_true(fit$converge)
  expect_equal(fit$fn, f_true, tolerance = 1e-4)
  expect_equal(dim(fit$B), c(n, 1L))

  # Single vector: up to sign flip
  dot <- abs(drop(crossprod(fit$B, v_true)))
  expect_equal(dot, 1, tolerance = 1e-4)
})

test_that("ortho_optim corrects non-orthogonal initial B", {
  set.seed(442)
  n <- 8
  k <- 2

  A <- matrix(rnorm(n * n), n, n)
  A <- crossprod(A)

  eig <- eigen(A, symmetric = TRUE)
  f_true <- -0.5 * sum(eig$values[1:k])

  # Deliberately non-orthogonal B
  B0 <- matrix(rnorm(n * k), n, k)

  fn <- function(B) -0.5 * sum(diag(crossprod(B, A %*% B)))
  gr <- function(B) -A %*% B

  # Gram-Schmidt correction is printed via cat(), not message()
  expect_output(
    fit <- ortho_optim(B0, fn, gr, maxitr = 500, verbose = FALSE),
    "Gram-Schmidt"
  )

  expect_true(fit$converge)
  expect_equal(fit$fn, f_true, tolerance = 1e-4)
})

test_that("ortho_optim validates inputs", {
  expect_error(ortho_optim(NULL, function(B) 0), "Initial value of B must be given")
  expect_error(ortho_optim(c(1, 2, 3), function(B) 0), "B must be a matrix")

  B <- diag(2)
  expect_error(ortho_optim(B, fn = NULL), "fn must be given")

  # fn must return single number
  expect_error(
    ortho_optim(B, fn = function(B) B),
    "fn must return a single number"
  )

  # grad must return matrix of same dimension
  expect_error(
    ortho_optim(B, fn = function(B) 1, grad = function(B) c(1, 2)),
    "grad must return a matrix"
  )
})

test_that("ortho_optim numerical and analytical gradients agree", {
  set.seed(542)
  n <- 8
  k <- 2

  A <- matrix(rnorm(n * n), n, n)
  A <- crossprod(A)

  B0 <- qr.Q(qr(matrix(rnorm(n * k), n, k)))

  fn <- function(B) -0.5 * sum(diag(crossprod(B, A %*% B)))
  gr <- function(B) -A %*% B

  fit_anal  <- ortho_optim(B0, fn, gr, maxitr = 500, verbose = FALSE)
  fit_numer <- ortho_optim(B0, fn, grad = NULL, maxitr = 500, verbose = FALSE)

  # Both should converge to same objective value
  expect_true(fit_anal$converge)
  expect_true(fit_numer$converge)
  expect_equal(fit_anal$fn, fit_numer$fn, tolerance = 1e-3)
})
