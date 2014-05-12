##
## function to find MSPE for cv data using cross-validated selection of sigma^2_beta using bayesian regression with shrinkage
##

fun.mcmc.chain.cv <- function(iter, data.cv, sigma.squared.beta){
  model.fit <- mcmc.lm.cv(data.cv[[iter]]$Y.cv, data.cv[[iter]]$X.cv, data.cv[[iter]]$Y.val, data.cv[[iter]]$X.val, n.mcmc, sigma.squared.beta)
  MSPE.cv <- mean((data.cv[[iter]]$Y.val - apply(model.fit$y.pred.save[, (n.burn + 1):n.mcmc], 1 , mean))^2)
  log.score.cv <-  mean(model.fit$log.score.save[(n.burn + 1):n.mcmc])
  return(c(MSPE.cv, log.score.cv))
}

make.grid.search.cv <- function(min.grid, max.grid, grid.size){
  sigma.squared.beta.cv <- seq(from = min.grid, to = max.grid, length = grid.size)
  MSPE.sim <- matrix(nrow = grid.size, ncol = 3)
  for(i in 1:grid.size){
    sigma.squared.beta <- sigma.squared.beta.cv[i]
    MSPE.sim[i, 1] <- sigma.squared.beta
    MSPE.sim[i, 2:3] <- apply(sfSapply(1:k.fold, fun.mcmc.chain.cv, data.cv = data.cv, sigma.squared.beta = sigma.squared.beta), 1, mean)
    cat(i, ' ')
  }
  return(MSPE.sim)
}