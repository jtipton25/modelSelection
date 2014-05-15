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

make.sim.data <- function(n, m, p, cor.vec, beta, sigma.squared.epsilon, sigma.squared.eta, phi){
  
  ## index of location
  locs.full <- 1:m
  
  ## predictor matrix of dimension p + 1 including intercept but not spatial location
  X.full <- matrix(nrow = m, ncol = p + 1)
  X.full[, 1:(p - 5 + 1)] <- matrix(c(rep(1, m), rnorm(m * (p - 5))), nrow = m, ncol = p - 5 + 1)
  X.full[, (p - 5 + 2):(p + 1)] <- X.full[, 2:6] %*% cor.vec %*% rep(1, 5) + matrix(rnorm(5 * n.o), nrow = m, ncol = 5)
  ## center and scale the predictor variables using the standard deviation denominator of m instead of m - 1
  scaled <- scale.predictor(X.full)
  X.full <- scaled$X
  scale <- scaled$scale 
  ## Construct spatial covariance matrix Sig.s
  D <- as.matrix(dist(1:m))
  Sig.s <- sigma.squared.eta * exp( - D / phi)
  ## Construct model covaraince matrix Sig
  Sig <- sigma.squared.epsilon * diag(m) + Sig.s
  
  ## simulate the response
  Y.full <- rmvnorm(1, X.full %*% beta, Sig)
  
  ##
  ## sample the spatial field
  ##
  
  samp <- sample(1:m, n)
  locs.o <- locs.full[samp]
  X.o <- X.full[samp, ]
  Y.o <- Y.full[samp]
    
  #   ##
  #   ## simulate new data for prediction
  #   ##
  # 
  #   X.new <- matrix(nrow = n.new, ncol = p + 1)
  #   X.new[, 1:(p - 5 + 1)] <- matrix(c(rep(1, n.new), rnorm(n.new * (p - 5))), nrow = n.new, ncol = p - 5 + 1)
  #   X.new[, (p - 5 + 2):(p + 1)] <- X.new[, 2:6] %*% cor.vec %*% rep(1, 5) + matrix(rnorm(5 * n.new), nrow = n.new, ncol = 5)
  # 
  #   ## center and scale the predictor variables using the standard deviation denomintor of n.o from above
  #   X.new.center <- matrix(nrow = n.new, ncol = p + 1) 
  #   X.new.center[, 1] <- rep(0, n.new)
  #   X.new.center[, 2:(p + 1)] <- (X.new[, 2:(p + 1)] - scale[,1]) / scale[, 2]
  #   X.new <- X.new.center  
  #   
  #   ## new Y values
  #   Y.new <- as.vector(rmvnorm(1, mean = X.new %*% beta, sigma = Sig.s))
  #   
  #   ## output the data
  #   list(X.o = X.o, Y.o = Y.o, X.new = X.new, Y.new = Y.new, scale = scale)
  # }
  ## output the data
  list(X.o = X.o, Y.o = Y.o, locs.o = locs.o, X.full = X.full, Y.full = Y.full, locs.full = locs.full, scale = scale)
} 