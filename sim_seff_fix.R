#!/usr/bin/env Rscript
# =============================================================================
# SEFF Fixes Verification: n=400, d=10, 4 models, 50 reps
# Fixes applied:
#   1. Y-kernel bandwidth: bw_y = N^{1/(ndr+5) - 1/7}  (was N^{-1/(ndr+6)})
#   2. initB: best-of-3 inits for SEFF/local (was SAVE-only)
# =============================================================================

.libPaths(c("~/R/library", .libPaths()))
library(orthoDr)

N      <- 400
P      <- 10
N_rep  <- 50
ncores <- 2
methods <- c("sir", "save", "local", "seff")

# True direction: 3 non-zero coefficients
B_true <- c(1, -0.5, 0.3, rep(0, P - 3))
B_true <- B_true / sqrt(sum(B_true^2))

# Heterosked: mean dir + variance dir
B_mean <- B_true
B_var  <- c(rep(0, 5), 1, -0.5, 0.3, rep(0, P - 8))
B_var  <- B_var / sqrt(sum(B_var^2))
B_het  <- cbind(B_mean, B_var)

settings <- list(
  list(label = "Linear",      ndr = 1, B = matrix(B_true, P, 1),
       gen = function(X, B) as.vector(X %*% B) + rnorm(N)),
  list(label = "Quadratic",   ndr = 1, B = matrix(B_true, P, 1),
       gen = function(X, B) as.vector((X %*% B)^2) + rnorm(N, 0, 0.5)),
  list(label = "Heterosked",  ndr = 2, B = B_het,
       gen = function(X, B) {
         u1 <- as.vector(X %*% B[,1]); u2 <- as.vector(X %*% B[,2])
         as.vector(rnorm(N, mean = u1, sd = 1 + abs(u2)))
       }),
  list(label = "Abs value",   ndr = 1, B = matrix(B_true, P, 1),
       gen = function(X, B) as.vector(abs(X %*% B)) + rnorm(N, 0, 0.5))
)

cat(sprintf("SEFF Fixes Verification\n"))
cat(sprintf("n=%d, d=%d, %d models x %d methods x %d reps, ncores=%d\n\n",
            N, P, length(settings), length(methods), N_rep, ncores))

results <- list()
t_all <- Sys.time()

for (si in seq_along(settings)) {
  s <- settings[[si]]
  cat(sprintf("=== %s (ndr=%d) ===\n", s$label, s$ndr))

  for (method in methods) {
    dists <- numeric(N_rep)
    times <- numeric(N_rep)
    convs <- integer(N_rep)
    errs  <- 0

    for (r in 1:N_rep) {
      X <- matrix(rnorm(N * P), N, P)
      Y <- s$gen(X, s$B)

      t0 <- Sys.time()
      tryCatch({
        fit <- orthoDr_reg(X, Y, method = method, ndr = s$ndr,
                           maxitr = 500, verbose = 0, ncore = ncores)
        times[r] <- as.numeric(Sys.time() - t0, units = "secs")
        dists[r] <- distance(s$B, fit$B, "dist")
        convs[r] <- fit$converge
      }, error = function(e) { errs <<- errs + 1 })
    }

    v <- times > 0
    if (sum(v) == 0) {
      cat(sprintf("  %-6s  ALL ERRORS\n", method))
      next
    }

    results[[paste(si, method, sep = "_")]] <- list(
      dist = dists[v], time = times[v], conv = convs[v], errs = errs
    )

    cat(sprintf("  %-6s  dist=%.4f±%.4f  time=%.3fs  conv=%d/%d%s\n",
                method, mean(dists[v]), sd(dists[v]), mean(times[v]),
                sum(convs[v]), sum(v),
                ifelse(errs > 0, sprintf(" err:%d", errs), "")))
  }
  cat("\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# Summary: Subspace Distance
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== SUBSPACE DISTANCE: Mean ± SD (lower=better) ===\n\n")

cat(sprintf("%-14s", "Model"))
for (m in methods) cat(sprintf("  %-14s", m))
cat("\n")
cat(strrep("-", 14 + 16 * length(methods)), "\n")

for (si in seq_along(settings)) {
  cat(sprintf("%-14s", settings[[si]]$label))
  for (method in methods) {
    r <- results[[paste(si, method, sep = "_")]]
    if (is.null(r)) { cat(sprintf("  %-14s", "FAIL")); next }
    cat(sprintf("  %6.4f±%5.4f", mean(r$dist), sd(r$dist)))
  }
  cat("\n")
}

# Convergence
cat("\n=== CONVERGENCE (out of 50) ===\n\n")
cat(sprintf("%-14s", "Model"))
for (m in methods) cat(sprintf("  %-14s", m))
cat("\n")
cat(strrep("-", 14 + 16 * length(methods)), "\n")

for (si in seq_along(settings)) {
  cat(sprintf("%-14s", settings[[si]]$label))
  for (method in methods) {
    r <- results[[paste(si, method, sep = "_")]]
    if (is.null(r)) { cat(sprintf("  %-14s", "FAIL")); next }
    cat(sprintf("  %-14s", sprintf("%d/%d", sum(r$conv), length(r$conv))))
  }
  cat("\n")
}

# Time
cat("\n=== MEAN TIME (s) ===\n\n")
cat(sprintf("%-14s", "Model"))
for (m in methods) cat(sprintf("  %-14s", m))
cat("\n")
cat(strrep("-", 14 + 16 * length(methods)), "\n")

for (si in seq_along(settings)) {
  cat(sprintf("%-14s", settings[[si]]$label))
  for (method in methods) {
    r <- results[[paste(si, method, sep = "_")]]
    if (is.null(r)) { cat(sprintf("  %-14s", "FAIL")); next }
    cat(sprintf("  %-14.2f", mean(r$time)))
  }
  cat("\n")
}

elapsed <- as.numeric(Sys.time() - t_all, "secs")
cat(sprintf("\nTotal: %.0fs (%.1f min)\n", elapsed, elapsed / 60))

saveRDS(results,
        file = "/home/tez/my-packages/orthoDr-project/orthoDr/sim_seff_fix_results.rds")
cat("Saved to sim_seff_fix_results.rds\n")
