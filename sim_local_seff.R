#!/usr/bin/env Rscript
# =============================================================================
# orthoDr Local & SEFF Benchmark
# =============================================================================
# Focused evaluation of local and seff methods across diverse settings.
# Also includes SIR and SAVE as reference methods.
#
# Settings cover the theoretical strengths/weaknesses of each method:
#   - Linear / monotone (SIR should win)
#   - Symmetric: quadratic, sine (local should fail, seff should handle)
#   - Asymmetric nonlinear: exp, |x|, mixture (local's strength)
#   - Heteroscedastic (local's other strength)
#   - Multi-index models with ndr=1 and ndr=2
# =============================================================================

library(orthoDr)
library(MASS)

N_rep  <- 50
ncores <- 2
methods <- c("sir", "save", "local", "seff")

# True B for ndr=1
B1 <- c(1, 1, 0, 0, 0, 0) / sqrt(2)

# True B for ndr=2
B2a <- c(1, 1, 0, 0, 0, 0) / sqrt(2)
B2b <- c(0, 0, 1, -1, 0, 0) / sqrt(2)
B2 <- cbind(B2a, B2b)

settings <- list(
  # ── ndr = 1, N = 200, P = 6 ──
  list(label = "S1: Linear",            ndr = 1, N = 200, P = 6, B = matrix(B1, 6, 1),
       gen = function(X, B) as.vector(X %*% B) + rnorm(nrow(X))),
  
  list(label = "S2: Quadratic",         ndr = 1, N = 200, P = 6, B = matrix(B1, 6, 1),
       gen = function(X, B) as.vector((X %*% B)^2) + rnorm(nrow(X), 0, 0.5)),
  
  list(label = "S3: Exponential",       ndr = 1, N = 200, P = 6, B = matrix(B1, 6, 1),
       gen = function(X, B) as.vector(exp(X %*% B)) + rnorm(nrow(X), 0, 0.5)),
  
  list(label = "S4: Abs value",         ndr = 1, N = 200, P = 6, B = matrix(B1, 6, 1),
       gen = function(X, B) as.vector(abs(X %*% B)) + rnorm(nrow(X), 0, 0.5)),
  
  list(label = "S5: Sine (symmetric)",  ndr = 1, N = 200, P = 6, B = matrix(B1, 6, 1),
       gen = function(X, B) as.vector(sin(pi * (X %*% B))) + rnorm(nrow(X), 0, 0.5)),
  
  list(label = "S6: Heterosc (local)",  ndr = 1, N = 200, P = 6, B = matrix(B1, 6, 1),
       gen = function(X, B) as.vector(rnorm(nrow(X), 0, 1 + (X %*% B)^2))),
  
  list(label = "S7: Cubic",             ndr = 1, N = 200, P = 6, B = matrix(B1, 6, 1),
       gen = function(X, B) as.vector((X %*% B)^3) + rnorm(nrow(X), 0, 0.5)),
  
  list(label = "S8: Mixture",           ndr = 1, N = 200, P = 6, B = matrix(B1, 6, 1),
       gen = function(X, B) {
         u <- as.vector(X %*% B)
         as.vector(ifelse(u > 0, u, -2*u)) + rnorm(nrow(X), 0, 0.5)
       }),
  
  # ── ndr = 1, N = 500 (higher N for harder settings) ──
  list(label = "S9:  Sine N=500",       ndr = 1, N = 500, P = 6, B = matrix(B1, 6, 1),
       gen = function(X, B) as.vector(sin(pi * (X %*% B))) + rnorm(nrow(X), 0, 0.5)),
  
  list(label = "S10: Heterosc N=500",   ndr = 1, N = 500, P = 6, B = matrix(B1, 6, 1),
       gen = function(X, B) as.vector(rnorm(nrow(X), 0, 1 + (X %*% B)^2))),
  
  # ── ndr = 2, N = 300, P = 6 ──
  list(label = "S11: Two linear (ndr2)", ndr = 2, N = 300, P = 6, B = B2,
       gen = function(X, B) as.vector(X %*% B[,1] + 0.5 * X %*% B[,2]) + rnorm(nrow(X))),
  
  list(label = "S12: Mean+Var (ndr2)",   ndr = 2, N = 300, P = 6, B = B2,
       gen = function(X, B) as.vector(X %*% B[,1] + (X %*% B[,2])^2) + rnorm(nrow(X), 0, 0.5)),
  
  list(label = "S13: Two quad (ndr2)",   ndr = 2, N = 300, P = 6, B = B2,
       gen = function(X, B) as.vector((X %*% B[,1])^2 + (X %*% B[,2])^2) + rnorm(nrow(X), 0, 0.5)),
  
  list(label = "S14: Lin+Exp (ndr2)",    ndr = 2, N = 300, P = 6, B = B2,
       gen = function(X, B) as.vector(X %*% B[,1] + exp(X %*% B[,2])) + rnorm(nrow(X), 0, 0.5))
)

cat(sprintf("orthoDr Local & SEFF Benchmark\n"))
cat(sprintf("%d settings × %d methods × %d reps, ncores=%d\n\n",
            length(settings), length(methods), N_rep, ncores))

results <- list()
t_all <- Sys.time()

for (si in seq_along(settings)) {
  s <- settings[[si]]
  cat(sprintf("=== %s ===\n", s$label))
  
  for (method in methods) {
    dists <- numeric(N_rep)
    times <- numeric(N_rep)
    convs <- integer(N_rep)
    errs  <- 0
    
    for (r in 1:N_rep) {
      X <- matrix(rnorm(s$N * s$P), s$N, s$P)
      Y <- s$gen(X, s$B)
      
      t0 <- Sys.time()
      tryCatch({
        fit <- orthoDr_reg(X, Y, method = method, ndr = s$ndr,
                           maxitr = 500, verbose = FALSE, ncore = ncores)
        times[r] <- as.numeric(Sys.time() - t0, units = "secs")
        dists[r] <- distance(s$B, fit$B, "dist")
        convs[r] <- fit$converge
      }, error = function(e) { errs <<- errs + 1 })
    }
    
    v <- times > 0
    if (sum(v) == 0) {
      cat(sprintf("  %-8s  ALL ERRORS\n", method))
      results[[paste(si, method, sep = "_")]] <- NULL
      next
    }
    
    results[[paste(si, method, sep = "_")]] <- list(
      dist = dists[v], time = times[v], conv = convs[v], errs = errs
    )
    
    cat(sprintf("  %-8s  dist=%.4f±%.4f  time=%.3fs  conv=%d/%d%s\n",
                method, mean(dists[v]), sd(dists[v]), mean(times[v]),
                sum(convs[v]), sum(v),
                ifelse(errs > 0, sprintf(" err:%d", errs), "")))
  }
  cat("\n")
}

# Summary table
cat("══════════════════════════════════════════════════════════════════════════════\n")
cat("     SUBSPACE DISTANCE: Mean ± SD   (lower = better)\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

cat(sprintf("%-28s", "Setting"))
for (m in methods) cat(sprintf("%-14s", m))
cat("\n")
cat(paste(rep("-", 28 + 14 * length(methods)), collapse = ""), "\n")

for (si in seq_along(settings)) {
  cat(sprintf("%-28s", settings[[si]]$label))
  for (method in methods) {
    r <- results[[paste(si, method, sep = "_")]]
    if (is.null(r)) { cat(sprintf("%-14s", "FAIL")); next }
    cat(sprintf("%.3f±%.3f    ", mean(r$dist), sd(r$dist)))
  }
  cat("\n")
}

# Convergence table
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("     CONVERGENCE RATE (out of 50)\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

cat(sprintf("%-28s", "Setting"))
for (m in methods) cat(sprintf("%-14s", m))
cat("\n")
cat(paste(rep("-", 28 + 14 * length(methods)), collapse = ""), "\n")

for (si in seq_along(settings)) {
  cat(sprintf("%-28s", settings[[si]]$label))
  for (method in methods) {
    r <- results[[paste(si, method, sep = "_")]]
    if (is.null(r)) { cat(sprintf("%-14s", "FAIL")); next }
    cat(sprintf("%d/%-10d", sum(r$conv), length(r$conv)))
  }
  cat("\n")
}

# Time table
cat("\n══════════════════════════════════════════════════════════════════════════════\n")
cat("     MEAN TIME (seconds)\n")
cat("══════════════════════════════════════════════════════════════════════════════\n\n")

cat(sprintf("%-28s", "Setting"))
for (m in methods) cat(sprintf("%-14s", m))
cat("\n")
cat(paste(rep("-", 28 + 14 * length(methods)), collapse = ""), "\n")

for (si in seq_along(settings)) {
  cat(sprintf("%-28s", settings[[si]]$label))
  for (method in methods) {
    r <- results[[paste(si, method, sep = "_")]]
    if (is.null(r)) { cat(sprintf("%-14s", "FAIL")); next }
    cat(sprintf("%-14.3f", mean(r$time)))
  }
  cat("\n")
}

cat(sprintf("\nTotal time: %.1f seconds (%.1f min)\n",
            as.numeric(Sys.time() - t_all, "secs"),
            as.numeric(Sys.time() - t_all, "mins")))

saveRDS(results,
        file = "/home/tez/my-packages/orthoDr-project/orthoDr/sim_local_seff_results.rds")
cat("Results saved to sim_local_seff_results.rds\n")
