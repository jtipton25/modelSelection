##
## Model selection using orthogonal data augmentation following Ghosh and Clyde: "Rao-blackwellization for Bayesian Variable Selection and Model Averaging in a Linear and Binary Regression: A Novel Data Augmentation Approach
##
set.seed(10)

## number of covariates
p <- 5

## number of observations
n.o <- 100

##
## libraries and functions
##

library(mvtnorm)


## predictor matrix of dimension p + 1 including intercept
X.o <- matrix(c(rep(1, n.o), rnorm(n.o * p)), nrow = n.o, ncol = p + 1)

## center and scale the predictor variables using the standard deviation denomintor of n.o
scale <- matrix(nrow = p, ncol = 2)
for(i in 1:p){
  scale[i, ] <- c(mean(X.o[, i + 1]), sqrt((n.o - 1) / n.o) * sd(X.o[, i + 1]))
  X.o[, i + 1] <- (X.o[, i + 1] - scale[i, 1]) / scale[i, 2]
}
  
## regression parameters
beta <- sample(1:20, 6)
sigma.squared <- 0.5

 Y.o <- X.o %*% beta + rnorm(n.o, mean = 0, sd = sqrt(sigma.squared))

##
## Data augmentation
##

## number of observations to augment
n.a
n.c <- n.o + n.c
