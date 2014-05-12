##
## function to find MSPE for cv data using cross-validated bayesian LASSO regression
##

fun.mcmc.chain.lasso <- function(iter, data.cv, lambda.squared, alpha.epsilon, beta.epsilon){
  model.fit <- mcmc.lm.lasso.cv(data.cv[[iter]]$Y.cv, data.cv[[iter]]$X.cv, data.cv[[iter]]$Y.val, data.cv[[iter]]$X.val, n.mcmc, alpha.epsilon, beta.epsilon, lambda.squared)
  MSPE.cv <- mean((data.cv[[iter]]$Y.val - apply(model.fit$y.pred.save[, (n.burn + 1):n.mcmc], 1 , mean))^2)
  log.score.cv <- mean(model.fit$log.score.save[(n.burn + 1):n.mcmc])
  return(c(MSPE.cv, log.score.cv))
}

make.grid.search.lasso <- function(min.grid, max.grid, grid.size){
  lambda.squared.cv <- seq(from = min.grid, to = max.grid, length = grid.size)
  MSPE.sim <- matrix(nrow = grid.size, ncol = 3)
  for(i in 1:grid.size){
    lambda.squared <- lambda.squared.cv[i]
    MSPE.sim[i, 1] <- lambda.squared
    MSPE.sim[i, 2:3] <- apply(sfSapply(1:k.fold, fun.mcmc.chain.lasso, data.cv = data.cv, lambda.squared = lambda.squared, alpha.epsilon = alpha.epsilon, beta.epsilon = beta.epsilon), 1, mean)
    cat(i, ' ')
  }
  return(MSPE.sim)
}
