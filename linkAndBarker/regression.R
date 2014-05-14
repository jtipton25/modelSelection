##
## Model selection using orthogonal data augmentation following Ghosh and Clyde: "Rao-blackwellization for Bayesian Variable Selection and Model Averaging in a Linear and Binary Regression: A Novel Data Augmentation Approach
##
rm(list = ls())
set.seed(10)

##
## libraries and subroutines
##

## simulate the data
source('~/modelSelection/functions/simNottKohn.R')
# source('~/modelSelection/functions/simulateData.R')'
## bayesian linear regression code
source('~/modelSelection/bayesianLinearRegression/mcmc.lm.R')
source('~/modelSelection/functions/dinvgamma.R')
source('~/modelSelection/functions/rMVN.R')
source('~/modelSelection/plots/make.model.plot.R')
##
## simulate the regression data
##

## number of covariates (not including the intercept)
p <- 15

## number of observations
n.o <- 100

## correlation vector used to correlate the first five and last five covariates\
cor.vec <- c(0.3, 0.5, 0.7, 0.9, 1.1)

## regression coefficients
beta <- c(4, 2, 0, 0, 0, -1, 0, 1.5, 0, 0, 0, 1, 0, 0.5, 0, 0)
sigma.squared <- 2.5
sigma.squared.epsilon <- sigma.squared

## number of new observations for out of sample validation
n.new <- 30

simdata <- make.sim.data(n.o, p, cor.vec, beta, sigma.squared, n.new)
X.o <- simdata$X.o
Y.o <- simdata$Y.o
X.new <- simdata$X.new
Y.new <- simdata$Y.new




# tau <- 8
# N <- 100
# n <- 30
# beta <- 1:tau
# sigma.squared.epsilon <- 1
# 
# data <- make.data(N, tau, beta, sigma.squared.epsilon)
# 
# samp <- sample(1:N, n)
# 
# data.samp <- data[samp, ]
# 
# summary(lm(y ~ ., data = data.samp))
# summary(lm(y ~ X1, data = data.samp))

##
## Setup priors
##

# hyperparameters for mu.beta and sigma.squared.beta
mu.0 <- c(0, 0)#rep(0, tau)
sigma.squared.beta.0 <- 100 
# hyerparamters for sigma.squared.beta
alpha.beta <- 2
beta.beta <- 10
curve(dinvgamma(x, alpha.beta, beta.beta), from = 0, to = 10)
# hyperparameters for sigma.squared.epsilon
alpha.epsilon <- 2
beta.epsilon <- 10
curve(dinvgamma(x, alpha.epsilon, beta.epsilon), from = 0, to = 10)
## prior model weights
pi <- rep(1 / p, p)
## mcmc iterations
n.mcmc <- 5000
n.burn <- n.mcmc / 5

##
## Fit mcmc
##


out <- vector('list', length = p)
for(i in 1:p){
  out[[i]] <- mcmc.lm(Y.o, cbind(X.o[, 1], X.o[, i + 1]), cbind(X.new[, 1], X.new[, i + 1]), n.mcmc, mu.0, sigma.squared.beta.0, alpha.beta, beta.beta, alpha.epsilon, beta.epsilon)
}

beta

make.model.plot(out[[1]])
make.model.plot(out[[5]])


make.posterior.model.weights <- function(out){
  n.burn <- n.mcmc / 5
  tau <- length(out)
  BIC <- vector(length = tau)
  for(i in 1:tau){
#     BIC[i] <- median(out[[i]]$BIC[(n.burn + 1):n.mcmc]) ## should we use posterior mean or median of BIC???
    BIC[i] <- mean(out[[i]]$BIC[(n.burn + 1):n.mcmc])
  }
  weights <- (exp( - BIC / 2) * pi) / sum(exp( - BIC / 2) * pi)
return(weights)
}

weights <- make.posterior.model.weights(out)
plot(weights)

tmp <- matrix(nrow = n.new, ncol = p)
for(i in 1:p){
  tmp[, i] <- apply(out[[i]]$y.pred.save[, (n.burn + 1):n.mcmc], 1, mean) * weights[i]
}
tmp
y.model.average <- apply(tmp, 1, sum)
MSPE.average <- mean((y.model.average - Y.new)^2)
MSPE.average


mu.0 <- rep(0, p + 1)
out.full <- mcmc.lm(Y.o, X.o, X.new, n.mcmc, mu.0, sigma.squared.beta.0, alpha.beta, beta.beta, alpha.epsilon, beta.epsilon)
MSPE.full <- mean((apply(out.full$y.pred.save, 1, mean) - Y.new)^2)

## seems that model averaging performs worse
MSPE.average
MSPE.full
