##
## MCMC sampler for orthogonal data augmentation
##

mcmc.oda <- function(Y.o, Y.a, X.o, X.a, params){
  
  ##
  ## functions and subroutines
  ##
  
  
  
  ##
  ## initialize values
  ##
 
  n.o <- length(Y.o)
  n.a <- length(Y.a)
  I.a <- diag(n.a)
  p <- dim(X.o)[2] - 1
  gamma <- sample(c(0, 1), p, replace = TRUE)
  X.o.gamma <- cbind(X.o[, 1], X.o[, 2:6][, gamma == 1])
  X.a.gamma <- cbind(X.a[, 1], X.a[, 2:6][, gamma == 1])
  alpha <- params[1]
  lambda <- c(0, rgamma(p, alpha / 2, alpha/ 2))
  Delta.gamma <- diag(c(0, lambda[2:6][which(gamma == 1)]))
  
  ##
  ## sample beta.tilde.gamma
  ##
  
  beta.tilde.gamma <- solve(t(X.o.gamma) %*% X.o.gamma + Delta.gamma) %*% t(X.o.gamma) %*% Y.o
  
  ##
  ## sample Y.a
  ##
  
  Y.a <- rmvnorm(1, mean = X.a.gamma %*% beta.tilde.gamma, sigma = sigma.squared * I.a + X.a.gamma %*% solve(t(X.o.gamma) %*% X.o.gamma + Delta.gamma) %*% t(X.a.gamma))
  
  ##
  ## sample sigma.squared
  ##
  
  sigma.squared <- 1 / rgamma(1, (n.o - 1) / 2, t(Y.o) %*% Y.o - t(beta.tilde.gamma) %*% (t(X.o.gamma) %*% X.o.gamma + Delta.gamma) %*% beta.tilde.gamma)
  
  
}