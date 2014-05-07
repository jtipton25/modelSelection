set.seed(10)

source('~/modelSelection/functions/simulateData.R')
source('~/modelSelection/bayesianLinearRegression/mcmc.lm.R')
source('~/modelSelection/functions/dinvgamma.R')
source('~/modelSelection/functions/rMVN.R')
source('~/modelSelection/plots/make.model.plot.R')

tau <- 8
N <- 100
n <- 30
beta <- 1:tau
sigma.squared.epsilon <- 1

data <- make.data(N, tau, beta, sigma.squared.epsilon)

samp <- sample(1:N, n)

data.samp <- data[samp, ]

summary(lm(y ~ ., data = data.samp))
summary(lm(y ~ X1, data = data.samp))

##
## Setup priors
##

# hyperparameters for mu.beta and sigma.squared.beta
mu.0 <- 0#rep(0, tau)
sigma.squared.beta.0 <- 100 
# hyerparamters for sigma.squared.beta
alpha.beta <- 2
beta.beta <- 10
curve(dinvgamma(x, alpha.beta, beta.beta), from = 0, to = 10)
# hyperparameters for sigma.squared.epsilon
alpha.epsilon <- 2
beta.epsilon <- 10
curve(dinvgamma(x, alpha.epsilon, beta.epsilon), from = 0, to = 10)
n.mcmc <- 5000

##
## Fit mcmc
##

Y <- data.samp[, 1]
X <- as.matrix(data.samp[, 2:(tau + 1)], ncol = tau)

out <- vector('list', length = tau)
for(i in 1:tau){
  out[[i]] <- mcmc.lm(Y, X[, i], data[, (i + 1)], n.mcmc, mu.0, sigma.squared.beta.0, alpha.beta, beta.beta, alpha.epsilon, beta.epsilon)
}

make.model.plot(out[[8]])

apply(out[[8]]$y.pred.save, 1, mean) - data[, 1 ]
y.mod.av <- vector(length = N)
for(i in 1:tau){
  y.mod.av <- y.mod.av + 1 / tau * apply(out[[i]]$y.pred.save, 1, mean)
}

y.mod.av - data[, 1]


mu.0 <- rep(0, 8)
out.full <- mcmc.lm(Y, X, as.matrix(data[, 2:9]), n.mcmc, mu.0, sigma.squared.beta.0, alpha.beta, beta.beta, alpha.epsilon, beta.epsilon)
apply(out.full$y.pred.save, 1, mean) - data[, 1]
