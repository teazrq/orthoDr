#!/usr/bin/env Rscript
library(orthoDr)

N_rep  <- 5
B1 <- c(1, 1, 0, 0, 0, 0) / sqrt(2)

settings <- list(
  list(label = "Linear",      gen = function(X, B) as.vector(X %*% B) + rnorm(nrow(X))),
  list(label = "Quadratic",   gen = function(X, B) as.vector((X %*% B)^2) + rnorm(nrow(X), 0, 0.5)),
  list(label = "Sine",        gen = function(X, B) as.vector(sin(pi * (X %*% B))) + rnorm(nrow(X), 0, 0.5)),
  list(label = "Exponential", gen = function(X, B) as.vector(exp(X %*% B)) + rnorm(nrow(X), 0, 0.5)),
  list(label = "Heterosc",    gen = function(X, B) as.vector(rnorm(nrow(X), 0, 1 + (X %*% B)^2))),
  list(label = "Abs value",   gen = function(X, B) as.vector(abs(X %*% B)) + rnorm(nrow(X), 0, 0.5))
)

methods <- c("sir", "save", "local", "seff")

for (s in settings) {
  cat(sprintf("\n=== %s (ndr=1, N=200, P=6) ===\n", s$label))
  for (method in methods) {
    dists <- numeric(N_rep)
    for (r in 1:N_rep) {
      X <- matrix(rnorm(200 * 6), 200, 6)
      Y <- s$gen(X, matrix(B1, 6, 1))
      tryCatch({
        fit <- orthoDr_reg(X, Y, method = method, ndr = 1, maxitr = 500, verbose = FALSE, ncore = 2)
        dists[r] <- distance(matrix(B1, 6, 1), fit$B, "sine")
      }, error = function(e) { dists[r] <- NA })
    }
    v <- !is.na(dists)
    cat(sprintf("  %-6s  sine = %.4f\n", method, mean(dists[v])))
  }
}
cat("\nDone.\n")
