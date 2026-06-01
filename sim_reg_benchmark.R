#!/usr/bin/env Rscript
# =============================================================================
# orthoDr Regression Benchmark: 6 settings × 5 methods × 50 reps
# =============================================================================

library(orthoDr)
library(MASS)

N_rep  <- 50
N      <- 150
P      <- 6
ncores <- 2
methods <- c("sir", "save", "phd", "seff", "local")

cat(sprintf("orthoDr Reg Benchmark: %d reps, N=%d, P=%d, ncores=%d\n",
            N_rep, N, P, ncores))

# True B matrices
B1 <- c(1, 1, 0, 0, 0, 0) / sqrt(2)
B2a <- c(1, 1, 0, 0, 0, 0) / sqrt(2)
B2b <- c(0, 0, 1, -1, 0, 0) / sqrt(2)
B2 <- cbind(B2a, B2b)

settings <- list(
  list(label = "S1: ndr=1, Linear",     ndr = 1, B = matrix(B1, P, 1),
       gen = function(X, B) 2*(X%*%B) + rnorm(nrow(X))),
  list(label = "S2: ndr=1, Quadratic",   ndr = 1, B = matrix(B1, P, 1),
       gen = function(X, B) (X%*%B)^2 + rnorm(nrow(X), 0, 0.5)),
  list(label = "S3: ndr=1, Exponential", ndr = 1, B = matrix(B1, P, 1),
       gen = function(X, B) exp(X%*%B) + rnorm(nrow(X), 0, 0.5)),
  list(label = "S4: ndr=2, Two linear",  ndr = 2, B = B2,
       gen = function(X, B) (X%*%B[,1]) + 0.5*(X%*%B[,2]) + rnorm(nrow(X))),
  list(label = "S5: ndr=2, Mean+Var",    ndr = 2, B = B2,
       gen = function(X, B) (X%*%B[,1]) + (X%*%B[,2])^2 + rnorm(nrow(X), 0, 0.5)),
  list(label = "S6: ndr=2, Two quad",    ndr = 2, B = B2,
       gen = function(X, B) (X%*%B[,1])^2 + (X%*%B[,2])^2 + rnorm(nrow(X), 0, 0.5))
)

sink("/home/tez/my-packages/orthoDr-project/orthoDr/sim_reg_output.txt")
cat(sprintf("orthoDr Regression Benchmark\n"))
cat(sprintf("N=%d, P=%d, reps=%d, ncores=%d\n\n", N, P, N_rep, ncores))

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
      X <- matrix(rnorm(N * P), N, P)
      Y <- as.vector(s$gen(X, s$B))
      
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
      results[[paste(si, method, sep="_")]] <- NULL
      next
    }
    
    results[[paste(si, method, sep="_")]] <- list(
      dist = dists[v], time = times[v], conv = convs[v], errs = errs
    )
    
    cat(sprintf("  %-8s  dist=%.4f±%.4f  time=%.3fs  conv=%d/%d%s\n",
                method, mean(dists[v]), sd(dists[v]), mean(times[v]),
                sum(convs[v]), sum(v),
                ifelse(errs > 0, sprintf(" err:%d", errs), "")))
  }
  cat("\n")
}

# Summary tables
cat("══════════════════════════════════════════════════════════════\n")
cat("     SUBSPACE DISTANCE: Mean ± SD   (lower = better)\n")
cat("══════════════════════════════════════════════════════════════\n\n")

cat(sprintf("%-28s", "Setting"))
for (m in methods) cat(sprintf("%-12s", m))
cat("\n")
cat(paste(rep("-", 28 + 12*length(methods)), collapse=""), "\n")

for (si in seq_along(settings)) {
  cat(sprintf("%-28s", settings[[si]]$label))
  for (method in methods) {
    r <- results[[paste(si, method, sep="_")]]
    if (is.null(r)) { cat(sprintf("%-12s", "FAIL")); next }
    cat(sprintf("%.3f±%.3f  ", mean(r$dist), sd(r$dist)))
  }
  cat("\n")
}

cat(sprintf("\nTotal time: %.1f seconds\n", as.numeric(Sys.time()-t_all, "secs")))
sink()

# Also print to stdout
cat(readLines("/home/tez/my-packages/orthoDr-project/orthoDr/sim_reg_output.txt"), sep="\n")

saveRDS(results, "/home/tez/my-packages/orthoDr-project/orthoDr/sim_reg_results.rds")
cat("\nResults saved.\n")
