##
## simulate regression model
##

make.data <- function(n, tau, beta, sigma.squared.epsilon){
  tau <- length(beta)
  X <- matrix(rnorm(n * tau), ncol = 8, nrow = n)
  y <- X %*% beta + rnorm(n, mean = 0, sd = sqrt(sigma.squared.epsilon))
  data.frame(y, X)  
}
