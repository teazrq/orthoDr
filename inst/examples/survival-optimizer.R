# Full survival and optimizer examples for revised development orthoDr 2.0.0.
# Uses FUN for plotting transformations; the CRAN 0.6.8 API also uses FUN.
# Source: inst/examples/survival-optimizer.R in teazrq/orthoDr.
# Run interactively for 3D rgl plots. This script is checked by R CMD check.

library(orthoDr)

# for the general solver
# an example of searching for the first principal component
library(pracma)
set.seed(1)
P = 100
N = 400

X = scale(matrix(rnorm(N*P), N, P), scale = FALSE)

w = matrix(rnorm(P), P, 1)
w = gramSchmidt(w)$Q

fx <- function(w, X) t(w) %*% t(X) %*% X %*% w
fit = ortho_optim(w, fx, X = X, maximize = TRUE, verbose = 0)
distance(fit$B, as.matrix(prcomp(X)$rotation[, 1]), type = "dist")

# if we have the gradient
gx <- function(w, X) 2*t(X) %*% X %*% w
fit = ortho_optim(w, fx, gx, X = X, maximize = TRUE, verbose = 2)
distance(fit$B, as.matrix(prcomp(X)$rotation[, 1]), type = "dist")



# for survival data
# Fitting the IR-CP model
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
ndr = 2

# fit the model
orthoDr.fit = orthoDr_surv(dataX, Y, Censor, ndr = ndr, verbose = 0)

# compare with the true direction
c(distance(failEDR, orthoDr.fit$B, "dist"),
   distance(failEDR, orthoDr.fit$B, "trace"),
   distance(failEDR, orthoDr.fit$B, "canonical", dataX))


# some plots
view_dr_surv(dataX, Y, Censor, orthoDr.fit$B[,1], FUN = log, type = "2D")
view_dr_surv(dataX, Y, Censor, orthoDr.fit$B[,1], FUN = log, type = "3D")


# the skcm data analysis with 20 genes

data(skcm.clinical)
data(skcm.melgene)

dataX = as.matrix(cbind(skcm.clinical$Age, skcm.melgene)[!is.na(skcm.clinical$Age), ])
dataX = scale(dataX)
Y = skcm.clinical$Time[!is.na(skcm.clinical$Age)]
Censor = skcm.clinical$Censor[!is.na(skcm.clinical$Age)]

orthoDr.fit = orthoDr_surv(dataX, Y, Censor, ndr = 1, verbose = 2)

CP_SIR(dataX, Y, Censor)$vectors[,1, drop = FALSE]

orthoDr.fit$B

view_dr_surv(dataX, sqrt(Y), Censor, orthoDr.fit$B[,1], FUN = identity, type = "2D")
view_dr_surv(dataX, sqrt(Y), Censor, -orthoDr.fit$B[,1], bw= 0.5, FUN = identity, type = "3D")
