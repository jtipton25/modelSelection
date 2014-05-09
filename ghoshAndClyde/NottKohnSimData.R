##
## Model selection using orthogonal data augmentation following Ghosh and Clyde: "Rao-blackwellization for Bayesian Variable Selection and Model Averaging in a Linear and Binary Regression: A Novel Data Augmentation Approach
##
rm(list = ls())
set.seed(10)

##
## libraries and subroutines
##

source('~/modelSelection/ghoshAndClyde/mcmc.ODA.R')
source('~/modelSelection/functions/rMVN.R')
library(snowfall)
library(rlecuyer)
library(parallel)

## number of covariates
p <- 15

## number of observations
n.o <- 50

## correlation vector used to correlate the first five and last five covariates\
cor.vec <- c(0.3, 0.5, 0.7, 0.9, 1.1)

##
## libraries and functions
##

library(mvtnorm)


## predictor matrix of dimension p + 1 including intercept
X.o <- matrix(nrow = n.o, ncol = p + 1)
X.o[, 1:(p - 5 + 1)] <- matrix(c(rep(1, n.o), rnorm(n.o * (p - 5))), nrow = n.o, ncol = p - 5 + 1)
X.o[, (p - 5 + 2):(p + 1)] <- X.o[, 2:6] %*% cor.vec %*% rep(1, 5) + matrix(rnorm(5 * n.o), nrow = n.o, ncol = 5)
  
## center and scale the predictor variables using the standard deviation denomintor of n.o
scale <- matrix(nrow = p, ncol = 2)
for(i in 1:p){
  scale[i, ] <- c(mean(X.o[, i + 1]), sqrt((n.o - 1) / n.o) * sd(X.o[, i + 1]))
  X.o[, i + 1] <- (X.o[, i + 1] - scale[i, 1]) / scale[i, 2]
}
  
## regression parameters
beta <- c(4, 2, 0, 0, 0, -1, 0, 1.5, 0, 0, 0, 1, 0, 0.5, 0, 0)
sigma.squared <- 2.5

 Y.o <- X.o %*% beta + rnorm(n.o, mean = 0, sd = sqrt(sigma.squared))

##
## Data augmentation
##

## number of observations to augment


alpha <- 2
pi.prior <- rep( 1 / 2, p)
epsilon = 0.001
n.mcmc <- 5000#200000
lambda <- c(0, rep(1, p))
params <- list('vector')
params <- list(n.mcmc, alpha, pi.prior, lambda)



## 
## fit mcmc
##

out <- mcmc.oda(Y.o = Y.o, X.o = X.o, params = params)


n.burn <- n.mcmc / 5
apply(out$beta.save, 1, mean)
beta.fit <- c(mean(out$beta.save[1, (n.burn + 1):n.mcmc]), apply(out$rho.save[, (n.burn + 1):n.mcmc] * out$delta.save / (out$delta.save + lambda[ - 1]) * out$beta.save[2:(p + 1), (n.burn + 1):n.mcmc], 1, mean))
    
      
    
    
##
## simulate new data for prediction
##

n.new <- 30
X.new <- matrix(nrow = n.new, ncol = p + 1)
X.new[, 1:(p - 5 + 1)] <- matrix(c(rep(1, n.new), rnorm(n.new * (p - 5))), nrow = n.new, ncol = p - 5 + 1)
X.new[, (p - 5 + 2):(p + 1)] <- X.new[, 2:6] %*% cor.vec %*% rep(1, 5) + matrix(rnorm(5 * n.new), nrow = n.new, ncol = 5)

## center and scale the predictor variables using the standard deviation denomintor of n.o
X.new.center <- matrix(nrow = n.new, ncol = p + 1) 
X.new.center[, 1] <- rep(0, n.new)
X.new.center[, 2:(p + 1)] <- (X.new[, 2:(p + 1)] - scale[,1]) / scale[, 2]

 ## regression parameters
beta <- c(4, 2, 0, 0, 0, -1, 0, 1.5, 0, 0, 0, 1, 0, 0.5, 0, 0)
sigma.squared <- 2.5

Y.new <- X.new.center %*% beta + rnorm(n.new, mean = 0, sd = sqrt(sigma.squared))
Y.new.hat <- X.new.center %*% beta.fit

## mean square prediction error
MSPE <- mean((Y.new - Y.new.hat)^2)

mod <- lm(Y.o ~ X.o[, 2:(p + 1)])
MSPE.lm <- mean((Y.new - X.new.center %*% coef(mod))^2)
MSPE
MSPE.lm

##
## Cross - validated selction on sigma^2_beta
## 


source('~/modelSelection/bayesianLinearRegressionCrossValidate/mcmc.lm.cv.R')
sigma.squared.beta <- 1
k.fold <- 8

make.cv.data <- function(Y, X, k.fold){
  n <- length(Y)
  p <- dim(X)[2]
  cv.samp.size <- rep(floor(n/ k.fold), k.fold) + c(rep(1, sum(n) %% k.fold), rep(0, k.fold - n %% k.fold))
  samp.cv <- vector('list', length = k.fold)
  out.cv <- vector('list', length = k.fold)
  for(iter in 1:k.fold){
    if(iter == 1){
      samp.cv[[iter]] <- sample(1:n, cv.samp.size[iter], replace = FALSE)  
    } else {
      tmp <- c()
      for(j in 1:iter){
        tmp <- c(tmp, samp.cv[[j]])
      }
    samp.cv[[iter]] <- sample((1:n)[ - tmp], cv.samp.size[iter], replace = FALSE)
    } 
    Y.cv <- Y[ - samp.cv[[iter]]]
    Y.val <- Y[samp.cv[[iter]]]
    X.cv <- X[ - samp.cv[[iter]], ]
    X.val <- X[samp.cv[[iter]], ]
    out.cv[[iter]] <- list(X.cv = X.cv, Y.cv = Y.cv, X.val = X.val, Y.val = Y.val)
  }
  return(out.cv)
}


data.cv <- make.cv.data(Y.o, X.o, k.fold)
#data.cv
n.burn <- n.mcmc / 5


## set up parallel cluster
sfInit(parallel = TRUE, cpus = 8)
sfExportAll()
sfClusterSetupRNG()
##
min.grid <- 0.25
max.grid <- 4
grid.size <- 16
sigma.squared.beta <- seq(from = min.grid, to = max.grid, length = grid.size)
sigma.squared.beta 

##
## functions
##
## function to find MSPE for cv data
fun.mcmc.chain <- function(iter, data.cv, sigma.squared.beta){
  model.fit <- mcmc.lm(data.cv[[iter]]$Y.cv, data.cv[[iter]]$X.cv, data.cv[[iter]]$X.val, n.mcmc, sigma.squared.beta)
  MSPE.cv <- mean((data.cv[[iter]]$Y.val - apply(model.fit$y.pred.save[, (n.burn + 1):n.mcmc], 1 , mean))^2)
  return(MSPE.cv)
}

make.grid.search <- function(min.grid, max.grid, grid.size){
  sigma.squared.beta <- seq(from = min.grid, to = max.grid, length = grid.size)
  MSPE.sim <- matrix(nrow = grid.size, ncol = 2)
  for(i in 1:grid.size){
    MSPE.sim[i, 1] <- sigma.squared.beta[i]
    MSPE.sim[i, 2] <- mean(sfSapply(1:k.fold, fun.mcmc.chain, data.cv = data.cv, sigma.squared.beta = sigma.squared.beta[i]))
    cat(i, ' ')
  }
  return(MSPE.sim)
}

test <- make.grid.search(min.grid, max.grid, grid.size)
plot(test[, 2] ~ test[, 1], type = 'l')



sfStop()

##
## fit simple model with cross-validated sigma.squared.beta
##
sigma.squared.beta <- test[, 1][which(test[, 2] == min(test[, 2]))]
out.cv <- mcmc.lm(Y.o, X.o, X.new.center, n.mcmc, sigma.squared.beta)
MSPE.cv <- mean((Y.new - apply(out.cv$y.pred.save[, (n.burn + 1):n.mcmc], 1 , mean))^2)

##
## fit using Lasso regression
##
source('~/modelSelection/bayesianLassoRegression/fixedEffectModel/mcmc.lm.lasso.fixed.lambda.R')

## function to find MSPE for cv data
fun.mcmc.chain.lasso <- function(iter, data.cv, lambda.squared, alpha.epsilon, beta.epsilon){
  model.fit <- mcmc.lm.lasso(data.cv[[iter]]$Y.cv, data.cv[[iter]]$X.cv, data.cv[[iter]]$X.val, n.mcmc, alpha.epsilon, beta.epsilon, lambda.squared)
  MSPE.cv <- mean((data.cv[[iter]]$Y.val - apply(model.fit$y.pred.save[, (n.burn + 1):n.mcmc], 1 , mean))^2)
  return(MSPE.cv)
}

make.grid.search.lasso <- function(min.grid, max.grid, grid.size){
  lambda.squared <- seq(from = min.grid, to = max.grid, length = grid.size)
  MSPE.sim <- matrix(nrow = grid.size, ncol = 2)
  for(i in 1:grid.size){
    MSPE.sim[i, 1] <- lambda.squared[i]
    MSPE.sim[i, 2] <- mean(sfSapply(1:k.fold, fun.mcmc.chain.lasso, data.cv = data.cv, lambda.squared = lambda.squared[i], alpha.epsilon = alpha.epsilon, beta.epsilon = beta.epsilon))
    cat(i, ' ')
  }
  return(MSPE.sim)
}

##
## Cross-validate Lasso for prediction
##
alpha.epsilon <- 0.01
beta.epsilon <- 0.01
alpha.lambda <- 1
beta.lambda <- 20

min.grid <- 0.25
max.grid <- 4
grid.size <- 16

sfInit(parallel = TRUE, cpus = 8)
sfExportAll()
sfClusterSetupRNG()
sfLibrary(statmod)

cross.validate.lasso <- make.grid.search.lasso(min.grid, max.grid, grid.size)


sfStop()


##
## fit simple model with cross-validated sigma.squared.beta
##

plot(cross.validate.lasso[, 2] ~ cross.validate.lasso[, 1], type = 'l')
lambda.squared <- cross.validate.lasso[, 1][which(cross.validate.lasso[, 2] == min(cross.validate.lasso[, 2]))]

out.lasso <- mcmc.lm.lasso(Y.o, X.o, X.new.center, n.mcmc, alpha.epsilon, beta.epsilon,lambda.squared)
MSPE.lasso <- mean((Y.new - apply(X.new.center %*% out.lasso$beta.save[, (n.burn + 1):n.mcmc], 1, mean))^2)
mean((Y.new - apply(out.lasso$y.pred.save[, (n.burn + 1):n.mcmc], 1, mean))^2)

MSPE.lasso


apply(X.new.center %*% out.lasso$beta.save[, (n.burn + 1):n.mcmc], 1, mean)
apply(out.lasso$y.pred.save[, (n.burn + 1):n.mcmc], 1, mean)


##
## compare MSPE from different methods
##
MSPE
MSPE.lm
MSPE.cv
MSPE.lasso

#save.image('~/modelSelection/data/ODAmcmc.RData')