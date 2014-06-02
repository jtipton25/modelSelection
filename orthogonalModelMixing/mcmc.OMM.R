##
## MCMC sampler for orthogonal data augmentation
##

mcmc.omm <- function(Y.o, X.o, Y.new, X.new, params, epsilon = 0.001){
  
  ##
  ## functions and subroutines
  ##
  
  
  
  ##
  ## initialize fixed values
  ##
  
  n.mcmc <- params[[1]]
  alpha <- params[[2]]
  pi.prior <- params[[3]]  
  lambda <- params[[4]]
  mu.0 <- params[[5]]
  v <- params[[6]]
  psi <- params[[7]]
  
  X.pca <- prcomp(X.o, center = TRUE, retx = TRUE)
  X <- X.pca$x
  lambda.pca <- X.pca$sdev^2
  Lambda <- diag(lambda.pca)
#   Lambda.determinant <- det(Lambda)
#   Lambda.inv <- solve(Lambda)
  n.o <- length(Y.o)
  p <- dim(X.o)[2]
  I.o <- diag(n.o)
  I.p <- diag(p)
  
  X.tilde <- rbind(X, I.p)
  Y.tilde <- c(Y.o, mu.0)

  d <- diag(t(X.tilde) %*% X.tilde)
  
  P.z <- vector('list', length = p)
  SSR.squared <- vector(length = p)
  for(i in 1:p){
    P.z[[i]] <- X.tilde[, i] %*% t(X.tilde[, i]) / as.numeric(t(X.tilde[, i]) %*% X.tilde[, i])
    SSR.squared[i] <- sum((P.z[[i]] %*% Y.tilde)^2)
  }
  
#   log.q.gamma <- sum(gamma * (log(pi.prior / (1 - pi.prior)) - 1 / 2 * log(d))) - (n.o + v) / 2 * log(v * psi + t(Y.tilde) %*% Y.tilde - sum(gamma * SSR.squared))

X.tilde.1 <- rep(1, n.o + 1)
Y.tilde.1 <- c(Y.o, mu.0[1])
SSR.squared.1 <- sum(((X.tilde.1 %*% t(X.tilde.1) / as.numeric(t(X.tilde.1) %*% X.tilde.1)) %*% Y.tilde.1)^2)

LN <- log(v * psi + t(Y.tilde) %*% Y.tilde - SSR.1)
LF <- log(v * psi + t(Y.tilde) %*% Y.tilde - sum(SSR.squared))

  num.trunc <- p
  tmp <- vector(length = num.trunc)
  for(j in 1:num.trunc){
    tmp[j] <- SSR.squared.1^j / (j * (v * psi + t(Y.tilde) %*% Y.tilde)^j)
  }
  LN.bar <- sum(tmp)
  LN.bar

  tmp <- vector(length = num.trunc)
  for(j in 1:num.trunc){
    tmp[j] <- sum(SSR.squared^j) / (j * (v * psi + t(Y.tilde) %*% Y.tilde)^j)
  }
  LF.bar <- sum(tmp)
  LF.bar



  tmp.pi <- matrix(nrow = p, ncol = num.trunc)
  for(j in 1:num.trunc){
    tmp.pi[, j] <- (n.o + v) * SSR.squared^j / (2 * j * (v * psi + t(Y.tilde) %*% Y.tilde)^j)
  }

p.i <- (pi.prior * exp( - 1 / 2 * log(d) + (LN - LF) / (LN.bar - LF.bar) * apply(tmp.pi, 1, sum))) / (1 - pi.prior + pi.prior * exp( - 1 / 2 * log(d) + (LN - LF) / (LN.bar - LF.bar) * apply(tmp.pi, 1, sum)))

for(i in 1:N){
  gamma <- rbinom(p, 1, p.i)
  mod <- lm(Y.o ~ X[, (1:p)[gamma == 1]] + 0)
  predict(mod, X.new)  
?predict

#    log.q.gamma <- sum(gamma * (log(pi.prior / (1 - pi.prior)) - 1 / 2 * log(d))) - (n.o + v) / 2 * log(v * psi + t(Y.tilde) %*% Y.tilde - sum(gamma * SSR.squared))


  ## initialize random values
  
  gamma <- rbinom(p, 1, pi.prior)
  #   lambda <- c(0, rgamma(p, alpha / 2, alpha/ 2))
#   lambda <- c(0, rep(1, p))
  Delta.gamma <- diag(c(0, lambda[2:(p + 1)][which(gamma == 1)]))
  
  ## subset based on gamma
  X.o.gamma <- cbind(X.o[, 1], X.o[, 2:(p + 1)][, gamma == 1])
  
  ## sample Y.a
  #  Y.a <- 
  # for first sample of Y.a
  
  
  ##
  ## orthogonal data augmentation
  ## 
  
  delta <- eigen(t(X.o) %*% X.o)$values[1]
  D <- (delta + epsilon) * I.a
  X.a <- chol(D - (t(X.o) %*% X.o))
  X.c <- rbind(X.o, X.a)
  projectXontoY <- solve(t(X.c) %*% X.c) %*% t(X.c)
  X.a.gamma <- cbind(X.a[, 1], X.a[, 2:(p + 1)][, gamma == 1])
  beta.tilde.gamma <- solve(t(X.o.gamma) %*% X.o.gamma + Delta.gamma) %*% t(X.o.gamma) %*% Y.o
  
  
  ##
  ## setup save variables
  ##
  
  gamma.save <- matrix(nrow = p, ncol = n.mcmc)
  sigma.squared.save <- vector(length = n.mcmc)
  beta.save <- matrix(nrow = p + 1, ncol = n.mcmc)
  rho.save <- matrix(nrow = p, ncol = n.mcmc)
  delta.save <- delta
  log.score.save <- vector(length = n.mcmc)  

  ##
  ## begin mcmc
  ##
  
  for(k in 1:n.mcmc){
    if(k %% 10000 == 0){
      cat(k, ' ')
    }
    
    
    ##
    ## sample Y.a
    ##
    
    Y.a <- rmvnorm(1, mean = X.a.gamma %*% beta.tilde.gamma, sigma = sigma.squared * (I.a + X.a.gamma %*% solve(t(X.o.gamma) %*% X.o.gamma + Delta.gamma) %*% t(X.a.gamma)))
    Y.c <- c(Y.o, Y.a)
    beta.hat <- projectXontoY %*% Y.c
    
    ##
    ## sample sigma.squared
    ##
    
    sigma.squared <- 1 / rgamma(1, (n.o - 1) / 2, t(Y.o) %*% Y.o - t(beta.tilde.gamma) %*% (t(X.o.gamma) %*% X.o.gamma + Delta.gamma) %*% beta.tilde.gamma)
    
    
    ##
    ## sample gammma
    ##
    
    O <- pi.prior / (1 - pi.prior) * (lambda[2:(p + 1)] / delta + lambda[2:(p + 1)])^(1 / 2) * exp(1 / 2 * delta / (delta + lambda[2:(p + 1)]) * beta.hat[2:(p + 1)]^2 / sigma.squared * delta) 
    rho <- O / (1 + O)
    gamma <- rbinom(p, 1, rho)
    Delta.gamma <- diag(c(0, lambda[2:(p + 1)][which(gamma == 1)]))
    
    ##
    ## sample beta.tilde.gamma
    ##
    
    X.o.gamma <- cbind(X.o[, 1], X.o[, 2:(p + 1)][, gamma == 1])
    X.a.gamma <- cbind(X.a[, 1], X.a[, 2:(p + 1)][, gamma == 1])
    beta.tilde.gamma <- solve(t(X.o.gamma) %*% X.o.gamma + Delta.gamma) %*% t(X.o.gamma) %*% Y.o
    
    ##
    ## log scoring rule
    ##
    
    log.score <- sum(dnorm(Y.new, mean = cbind(X.new[, 1], X.new[, 2:(p + 1)][, gamma == 1]) %*% beta.tilde.gamma, sd = sqrt(sigma.squared), log = TRUE))
    
    ##
    ## save samples
    ##
    
    gamma.save[, k] <- gamma
    sigma.squared.save[k] <- sigma.squared
    beta.save[, k] <- beta.hat 
    rho.save[, k] <- rho
    delta.save <- delta
    log.score.save[k] <- log.score
  }
  list(gamma.save = gamma.save, sigma.squared.save = sigma.squared.save, beta.save = beta.save, rho.save = rho.save, delta.save = delta.save, log.score.save = log.score.save)
  
}

