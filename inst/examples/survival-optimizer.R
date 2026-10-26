# Supplemental examples for orthoDr 2.0.0: principal components, a simulated
# survival model, and the bundled melanoma data. Run interactively to explore
# the 3D rgl plots. The simpler survival tutorial plots predicted curves.

library(orthoDr)

# Find the first principal component using the general optimizer.
library(pracma)
set.seed(1)
P = 100
N = 400

X = scale(matrix(rnorm(N*P), N, P), scale = FALSE)

w = matrix(rnorm(P), P, 1)
w = gramSchmidt(w)$Q

fx <- function(w, X) t(w) %*% t(X) %*% X %*% w
fit = ortho_optim(w, fx, X = X, maximize = TRUE, verbose = 0)
c(subspace_distance = distance(fit$B, as.matrix(prcomp(X)$rotation[, 1]), type = "dist"))
# A distance near zero indicates agreement with the principal component.

# Providing the analytic gradient solves the same optimization problem.
gx <- function(w, X) 2*t(X) %*% X %*% w
fit = ortho_optim(w, fx, gx, X = X, maximize = TRUE, verbose = 0)
c(subspace_distance = distance(fit$B, as.matrix(prcomp(X)$rotation[, 1]), type = "dist"))



# Fit a survival model with two reduced predictors.
# In this simulation, event times depend on two linear combinations of X.
set.seed(1)
P = 6
N = 350
dataX = matrix(rnorm(N*P), N, P)
failEDR = as.matrix(cbind(c(1, 1, 0, 0, 0, 0, rep(0, P-6)),
                          c(0, 0, 1, -1, 0, 0, rep(0, P-6))))
censorEDR = as.matrix(c(0, 1, 0, 1, 1, 1, rep(0, P-6)))
T = exp( -2.5 + dataX %*% failEDR[,1] + 0.5 * (dataX %*% failEDR[,1]) * (dataX %*% failEDR[,2]) + 0.25 * log(-log(1-runif(N)))  )
C = exp( -0.5 + dataX %*% censorEDR  + log(-log(1-runif(N))) )
Y = pmin(T, C)
Censor = (T < C)
# Despite the argument name, TRUE (1) indicates an event; FALSE (0) is censoring.
ndr = 2

# Use the martingale estimating equation and request the two simulated directions.
orthoDr.fit = orthoDr_surv(dataX, Y, Censor, method = "dm", ndr = ndr, verbose = 0)

# Compare the estimated subspace with the known simulated subspace.
# Distance is better near 0; both correlations are better near 1.
c(subspace_distance = distance(failEDR, orthoDr.fit$B, "dist"),
  trace_correlation = distance(failEDR, orthoDr.fit$B, "trace"),
  canonical_correlation = distance(failEDR, orthoDr.fit$B, "canonical", dataX))


# Each plot below shows only the first coordinate of the two-dimensional fit.
# It does not display the full fitted subspace. FUN = log transforms time.
view_dr_surv(dataX, Y, Censor, orthoDr.fit$B[,1], FUN = log, type = "2D",
             ylab = "Log observed time")
view_dr_surv(dataX, Y, Censor, orthoDr.fit$B[,1], FUN = log, type = "3D",
             ylab = "Log observed time")


# Summarize survival in the bundled TCGA cutaneous melanoma (SKCM) data.
# Use age and 20 gene-expression predictors, excluding observations with
# missing age, and standardize these predictor columns before fitting.
# Here we request one direction to obtain a one-dimensional summary; this is
# an analyst choice, not an estimated structural dimension.

data(skcm.clinical)
data(skcm.melgene)

dataX = as.matrix(cbind(skcm.clinical$Age, skcm.melgene)[!is.na(skcm.clinical$Age), ])
dataX = scale(dataX)
Y = skcm.clinical$Time[!is.na(skcm.clinical$Age)]
Censor = skcm.clinical$Censor[!is.na(skcm.clinical$Age)]

orthoDr.fit = orthoDr_surv(dataX, Y, Censor, ndr = 1, verbose = 0)

# Compare with a CP-SIR estimate from the same observations.
# The distance compares subspaces without depending on the direction's sign.
cp_direction = CP_SIR(dataX, Y, Censor)$vectors[,1, drop = FALSE]
c(subspace_distance = distance(cp_direction, orthoDr.fit$B, "dist"))

orthoDr.fit$B

# Time is already square-root transformed, so FUN = identity leaves it as is.
# The sign flip in the 3D view changes the display orientation, not the subspace.
view_dr_surv(dataX, sqrt(Y), Censor, orthoDr.fit$B[,1], FUN = identity, type = "2D",
             ylab = "Square root of observed time")
view_dr_surv(dataX, sqrt(Y), Censor, -orthoDr.fit$B[,1], bw= 0.5, FUN = identity, type = "3D",
             ylab = "Square root of observed time")
