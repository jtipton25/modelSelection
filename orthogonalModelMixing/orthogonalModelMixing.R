##
## Model selection using orthogonal data augmentation following Ghosh and Clyde: "Rao-blackwellization for Bayesian Variable Selection and Model Averaging in a Linear and Binary Regression: A Novel Data Augmentation Approach
##
rm(list = ls())
set.seed(10)

##
## libraries and subroutines
##

source('~/modelSelection/functions/rMVN.R')
## simulate the data
source('~/modelSelection/functions/simNottKohn.R')
## make the cross-validation data
# source('~/modelSelection/functions/make.cv.data.R')
## load the OMM mcmc code
source('~/modelSelection/orthogonalModelMixing/mcmc.OMM.R')
## load shrinkage mcmc code
# source('~/modelSelection/bayesianLinearRegressionCrossValidate/mcmc.lm.cv.R') # shrinkage model with shrinkage parameter chosen by cross-validation
# source('~/modelSelection/functions/make.shrinkage.cv.R') # wrapper for parallelization of cross-validation of shrinkage model
## load lasso mcmc code
# source('~/modelSelection/bayesianLassoRegression/fixedEffectModel/mcmc.lm.lasso.R') # lasso with shrinkage parameter in model
# source('~/modelSelection/bayesianLassoRegression/fixedEffectModel/mcmc.lm.lasso.fixed.lambda.R') # lasso with shrinkage parameter fixed by cross-validation
# source('~/modelSelection/functions/make.lasso.cv.R') # wrapper for paralleliztion of cross-validation of lasso model

library(statmod)
library(mvtnorm)

##
## simulate the regression data
##


## number of covariates (including the intercept)
p <- 16

## number of observations
n.o <- 50

## correlation vector used to correlate the first five and last five covariates\
cor.vec <- c(0.3, 0.5, 0.7, 0.9, 1.1)

## regression coefficients
beta <- c(4, 2, 0, 0, 0, -1, 0, 1.5, 0, 0, 0, 1, 0, 0.5, 0, 0)
sigma.squared <- 2.5

## number of new observations for out of sample validation
n.new <- 30

simdata <- make.sim.data(n.o, p, cor.vec, beta, sigma.squared, n.new)
X.o <- simdata$X.o
Y.o <- simdata$Y.o
X.new <- simdata$X.new
Y.new <- simdata$Y.new

##
## Data augmentation
##

alpha <- 2
pi.prior <- rep( 1 / 2, p)
epsilon = 0.001
n.mcmc <- 5000
lambda <- c(0, rep(1, p))
mu.0 <- rep(0, p)
v <- 5
psi <- 3
params <- list('vector')
params <- list(n.mcmc, alpha, pi.prior, lambda, mu.0, v, psi)
n.burn <- n.mcmc / 5
k.fold <- 8

## 
## fit mcmc using ODA model
##

out <- mcmc.omm(Y.o = Y.o, X.o = X.o, Y.new, X.new, params = params)

## Rao-blackwell estimates
beta.fit <- c(mean(out$beta.save[1, (n.burn + 1):n.mcmc]), apply(out$rho.save[, (n.burn + 1):n.mcmc] * out$delta.save / (out$delta.save + lambda[ - 1]) * out$beta.save[2:(p + 1), (n.burn + 1):n.mcmc], 1, mean))
Y.new.hat <- X.new %*% beta.fit

## mean square prediction error
MSPE <- mean((Y.new - Y.new.hat)^2)
log.score <- mean(out$log.score.save[(n.burn + 1):n.mcmc])
