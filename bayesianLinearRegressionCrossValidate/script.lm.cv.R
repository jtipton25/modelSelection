##
## libraries and functions
##

source('~/functions/dinvgamma.R')
source('~/Linear-Model/Simple Linear Regression/mcmc.lm.R')
source('~/functions/rMVN.R')

setwd("~/Linear-Model/Simple Linear Regression/")
source("mcmc.lm.R")

make.model.plot <- function(out){
  n.burn <- floor(n.mcmc / 5) + 1
  layout(matrix(1:4, 2))
  matplot(t(out$beta.save[, n.burn:n.mcmc]), type = 'l')
  #hist(out$beta.save[1,][n.burn:n.mcmc])
  #abline(v = beta[1], col = 'red')
  hist(out$beta.save[2,][n.burn:n.mcmc])
  abline(v = beta[2], col = 'red')
  plot(out$sigma.squared.beta.save[n.burn:n.mcmc], type = 'l')
  plot(out$sigma.squared.epsilon.save[n.burn:n.mcmc], type = 'l')
  abline(h = sigma.squared.epsilon, col = 'red')
}


##
## Simulate some data
##

N <- 1000
n <- 100
beta <- -3:3
sigma.squared.epsilon <- 0.25
tau <- length(beta)

make.lm.data <- function(N, n, beta, sigma.sqaured.epsilon){
  tau <- length(beta)
  X <- matrix(nrow = N, ncol = tau)
  for(i in 1:tau){
    X[, i] <- rnorm(N, 0, 1)
  }
  #   X <-matrix(c(rep(1, N), rep(seq(0, 1, length = N), tau - 1)), nrow = N, ncol = tau)
  #   if(is.null(dim(X))){
  #     Y <- X * beta + rnorm(N, 0, sigma.squared.epsilon)
  #   } else {
  Y <- X %*% beta + rnorm(N, 0, sigma.squared.epsilon)
  #   }
  #list(X = X, Y = Y, N = N, n = n, sigma.squared.epsilon = sigma.squared.epsilon)
  data.frame(Y, X)
}

data <- make.lm.data(N, n, beta, sigma.squared.epsilon)

samp <- sample(1:N, n)
data.samp <- data[samp, ]

lm(Y ~ . ,data = data)

##
## Setup priors
##

# hyperparameters for mu.beta and sigma.squared.beta
mu.0 <- rep(0, tau)
sigma.squared.0 <- 100 
# hyerparamters for sigma.squared.beta
alpha.beta <- 2
beta.beta <- 10
curve(dinvgamma(x, alpha.beta, beta.beta), from = 0, to = 10)
# hyperparameters for sigma.squared.epsilon
alpha.epsilon <- 2
beta.epsilon <- 10
curve(dinvgamma(x, alpha.epsilon, beta.epsilon), from = 0, to = 10)
n.mcmc <- 5000

##
## Fit mcmc
##

Y <- data.samp[, 1]
X <- as.matrix(data.samp[, 2:(tau + 1)], ncol = tau)

out <- mcmc.lm(Y, X, n.mcmc, mu.0, sigma.squared.beta.0, alpha.beta, beta.beta, alpha.epsilon, beta.epsilon)

make.model.plot(out)

