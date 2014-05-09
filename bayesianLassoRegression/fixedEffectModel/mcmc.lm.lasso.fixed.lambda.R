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

mcmc.lm.lasso <- function(Y, X, X.tilde, n.mcmc, alpha.epsilon, beta.epsilon, lambda.squared){
  
  ##
  ## Initialize variables
  ##
  
  n <- length(Y)
  n.val <- dim(X.tilde)[1]
  tau <- dim(X)[2]
  # sigma.squared.epsilon
  sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon, beta.epsilon) 
  # lambda
  #   lambda.squared <- rgamma(1, alpha.lambda, beta.lambda)
  # gamma
  gamma.squared <- rgamma(tau, 1, lambda.squared / 2)
  Dgamma <- diag(gamma.squared)
  # beta
  mu.beta <- rep(0, tau)
  beta <- rMVN(chol(sigma.squared.epsilon * Dgamma), mu.beta)
  
  ##
  ## save variables
  ##
  
  beta.save <- matrix(nrow = tau, ncol = n.mcmc)
  sigma.squared.epsilon.save <- vector(length = n.mcmc)
  gamma.squared.save <-  matrix(nrow = tau, ncol = n.mcmc)
  #   lambda.squared.save <- vector(length = n.mcmc)
  Dbar.save <- vector(length = n.mcmc)
  y.pred.save <-  matrix(nrow = n.val, ncol = n.mcmc)
  
  ##
  ## Start MCMC
  ##
  
  for(k in 1:n.mcmc){
    if(k %% 100 == 0){
      cat(" ", k)
    }
    
    ##
    ## sample beta
    ##
    
    A.chol <- chol(t(X) %*% X + Dgamma)
    b <- (t(X) %*% Y)
    beta <- rMVN(A.chol, b)
    
    ##
    ## sample sigma.squared.epsilon
    ##
    
    sigma.squared.epsilon <- 1 / rgamma(1, alpha.epsilon + n / 2 + tau / 2, beta.epsilon + 1 / 2 * t(Y - X %*% beta) %*% (Y - X %*% beta) + 1 / 2 * t(beta) %*% Dgamma %*% beta)
    
    ##
    ## sample gamma.squared
    ##
    
    mu.tilde <- sqrt(lambda.squared * sigma.squared.epsilon / beta^2)
    lambda.tilde <- lambda.squared
    gamma.squared <- rinvgauss(tau, mu.tilde, lambda.tilde)
    Dgamma <- diag(gamma.squared[,1])
    
    ##
    ## lambda.squared
    ##
    
    #     lambda.squared <- rgamma(1, alpha.lambda + tau, beta.lambda + sum(gamma.squared) / 2)
    
    ##
    ## DIC calculations
    ##
    
    Dbar.save[k] <- - 2 * sum(dnorm(Y, X %*% beta, sqrt(sigma.squared.epsilon), log = TRUE))
    
    ##
    ## Posterior Predictive
    ##
    
    y.pred <- X.tilde %*% beta + rnorm(n.val, mean = 0, sd = sqrt(sigma.squared.epsilon))
    
    ##
    ## save variables
    ##
    
    beta.save[, k] <- beta
    sigma.squared.epsilon.save[k] <- sigma.squared.epsilon
    gamma.squared.save[, k] <- gamma.squared
    #     lambda.squared.save[k] <- lambda.squared
    y.pred.save[, k] <- y.pred 
  }
  
  ##
  ## output
  ##
  
  list(beta.save = beta.save, sigma.squared.epsilon.save = sigma.squared.epsilon.save, gamma.squared.save = gamma.squared.save, y.pred.save = y.pred.save)
}