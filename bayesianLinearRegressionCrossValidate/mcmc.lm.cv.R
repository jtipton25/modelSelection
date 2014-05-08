##
## Simple linear regression model 
##
## John Tipton - created 01.25.2014
##

##
## model: Y = X %*% beta + epsilon
##

##
## libraries and functions
##

mcmc.lm <- function(Y, X, X.tilde, n.mcmc, sigma.squared.beta){
  
  ##
  ## Initialize variables
  ##
  
  n <- length(Y)
  if(is.null(dim(X))){tau <- 1} else {tau <- dim(X)[2]}
  if(is.null(dim(X.tilde))){N <- length(X.tilde)} else {N <- dim(X.tilde)[1]}
  #
  I.beta <- diag(tau)
  #   sigma.squared.beta <- 1 / rgamma(1, alpha.beta, beta.beta)
  Sigma.beta <- sigma.squared.beta * I.beta
  Sigma.beta.inv <- 1 / sigma.squared.beta * I.beta 
  #
  I.epsilon <- diag(n)
  I.full <- diag(N)
  #
  
  
## fixed effects
  beta <- rMVN(A.chol = chol(Sigma.beta), b = rep(0, tau))

## random effects
#   Sigma.0 <- sigma.squared.beta.0 * I.beta
#   Sigma.0.inv <- 1 / sigma.squared.beta.0 * I.beta
#   mu.beta <- rMVN(A.chol = chol(Sigma.0), b = Sigma.0 %*% mu.0)
#   beta <- rMVN(A.chol = chol(Sigma.beta), b = Sigma.beta %*% mu.beta) 
  
  ##
  ## save variables
  ##
  
  beta.save <- matrix(nrow = tau, ncol = n.mcmc)
#   sigma.squared.beta.save <- vector(length = n.mcmc)
  sigma.squared.epsilon.save <- vector(length = n.mcmc)
  Dbar.save <- vector(length = n.mcmc)
  y.pred.save <- matrix(nrow = N, ncol = n.mcmc)
  
  ##
  ## Start MCMC
  ##
  
  for(k in 1:n.mcmc){
    if(k %% 100 == 0){
      cat(" ", k)
    }
    
    
    ##
    ## sample sigma.squared.epsilon
    ##
    
    sigma.squared.epsilon <- 1 / rgamma(1, n / 2, 1 / 2 * t(Y - X %*% beta) %*% (Y - X %*% beta))
    Sigma.epsilon <- sigma.squared.epsilon * I.epsilon
    Sigma.epsilon.inv <- 1 / sigma.squared.epsilon* I.epsilon
    
    
    ##
    ## sample beta - fixed effects
    ##
        
        A.chol <- chol(t(X) %*% Sigma.epsilon.inv %*% X + Sigma.beta.inv)
        b <- (t(X) %*% Sigma.epsilon.inv %*% Y)
        beta <- rMVN(A.chol, b)
    
    ##
    ## sample beta - random effects
    ##
    #     
    #     A.chol <- chol(t(X) %*% Sigma.epsilon.inv %*% X + Sigma.beta.inv)
    #     b <- (t(X) %*% Sigma.epsilon.inv %*% Y + Sigma.beta.inv %*% mu.beta)
    #     beta <- rMVN(A.chol, b)
    
    ##
    ## sample mu.beta - random effects
    ##
    #     
    #     A.chol <- chol(Sigma.beta.inv + Sigma.0.inv)
    #     b <- (Sigma.beta.inv %*% beta + Sigma.0.inv %*% mu.0)
    #     mu.beta <- rMVN(A.chol, b)
    
    ##
    ## sample sigma.squared.beta - fixed for cross validation
    ##
    #     
    #     sigma.squared.beta <- 1 / rgamma(1, alpha.beta + tau / 2, beta.beta + 1 / 2 * t(beta - mu.beta) %*% (beta - mu.beta))
    #     Sigma.beta <- sigma.squared.beta * I.beta
    #     Sigma.beta.inv <- 1 / sigma.squared.beta * I.beta
    

    
    ##
    ## DIC calculations
    ##
    
    Dbar.save[k] <- - 2 * sum(dnorm(Y, X %*% beta, sqrt(sigma.squared.epsilon), log = TRUE))
    
    ##
    ## posterior predictive distribution
    ##
    
    y.pred.save[, k] <- rMVN(1 / sqrt(sigma.squared.epsilon) * I.full, 1 / sigma.squared.epsilon * X.tilde %*% beta)
    
    ##
    ## save variables
    ##
    
    beta.save[, k] <- beta
    sigma.squared.epsilon.save[k] <- sigma.squared.epsilon
#     sigma.squared.beta.save[k] <- sigma.squared.beta
  }
  
  ##
  ## output
  ##
  
  list(beta.save = beta.save, sigma.squared.epsilon.save = sigma.squared.epsilon.save, Dbar.save = Dbar.save, y.pred.save = y.pred.save)
}