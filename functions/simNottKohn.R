##
## Simulate data using Nott and Kohn simulation
##

##
## functions and subroutines
##

##  scale predictors using standard deviation divided by n instead of n-1

scale.predictor <- function(X){
  n <- dim(X)[1]
  p <- dim(X)[2] - 1
  scale <- matrix(nrow = p, ncol = 2)
  X.tmp <- X
  for(i in 1:p){
    scale[i, ] <- c(mean(X[, i + 1]), sqrt((n - 1) / n) * sd(X[, i + 1]))
    X.tmp[, i + 1] <- (X[, i + 1] - scale[i, 1]) / scale[i, 2]
  }
  list(X = X.tmp, scale = scale)
}

##
## function to simulate data
##

make.sim.data <- function(n.o, p, cor.vec, beta, sigma.squared, n.new){

  ## predictor matrix of dimension p + 1 including intercept
  X.o <- matrix(nrow = n.o, ncol = p + 1)
  X.o[, 1:(p - 5 + 1)] <- matrix(c(rep(1, n.o), rnorm(n.o * (p - 5))), nrow = n.o, ncol = p - 5 + 1)
  X.o[, (p - 5 + 2):(p + 1)] <- X.o[, 2:6] %*% cor.vec %*% rep(1, 5) + matrix(rnorm(5 * n.o), nrow = n.o, ncol = 5)

  scaled <- scale.predictor(X.o)
  X.o <- scaled$X
  scale <- scaled$scale
  ## center and scale the predictor variables using the standard deviation denomintor of n.o

  Y.o <- X.o %*% beta + rnorm(n.o, mean = 0, sd = sqrt(sigma.squared))

  ##
  ## simulate new data for prediction
  ##

  X.new <- matrix(nrow = n.new, ncol = p + 1)
  X.new[, 1:(p - 5 + 1)] <- matrix(c(rep(1, n.new), rnorm(n.new * (p - 5))), nrow = n.new, ncol = p - 5 + 1)
  X.new[, (p - 5 + 2):(p + 1)] <- X.new[, 2:6] %*% cor.vec %*% rep(1, 5) + matrix(rnorm(5 * n.new), nrow = n.new, ncol = 5)

  ## center and scale the predictor variables using the standard deviation denomintor of n.o from above
  X.new.center <- matrix(nrow = n.new, ncol = p + 1) 
  X.new.center[, 1] <- rep(0, n.new)
  X.new.center[, 2:(p + 1)] <- (X.new[, 2:(p + 1)] - scale[,1]) / scale[, 2]
  X.new <- X.new.center  
  
  ## new Y values
  Y.new <- X.new.center %*% beta + rnorm(n.new, mean = 0, sd = sqrt(sigma.squared))
  
  ## output the data
  list(X.o = X.o, Y.o = Y.o, X.new = X.new, Y.new = Y.new, scale = scale)
}