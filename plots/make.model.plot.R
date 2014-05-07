make.model.plot <- function(out){
  n.burn <- floor(n.mcmc / 5) + 1
  layout(matrix(1:4, 2))
  if(dim(out$beta.save)[1] == 1){
    matplot(out$beta.save[n.burn:n.mcmc], type = 'l')
  } else {
    matplot(t(out$beta.save[, n.burn:n.mcmc]), type = 'l')
  }
  #hist(out$beta.save[1,][n.burn:n.mcmc])
  #abline(v = beta[1], col = 'red')
#   hist(out$beta.save[2,][n.burn:n.mcmc])
#   abline(v = beta[2], col = 'red')
  plot(out$sigma.squared.beta.save[n.burn:n.mcmc], type = 'l')
  plot(out$sigma.squared.epsilon.save[n.burn:n.mcmc], type = 'l')
  abline(h = sigma.squared.epsilon, col = 'red')
}