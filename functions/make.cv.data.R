##
## function to make data for cross-validation
## 
## takes vector Y and matrix X and makes k-fold fit and validation data sets


make.cv.data <- function(Y, X, k.fold){
  n <- length(Y)
  p <- dim(X)[2]
  cv.samp.size <- rep(floor(n / k.fold), k.fold) + c(rep(1, sum(n) %% k.fold), rep(0, k.fold - n %% k.fold))
  samp.cv <- vector('list', length = k.fold)
  out.cv <- vector('list', length = k.fold)
  for(iter in 1:k.fold){
    if(iter == 1){
      samp.cv[[iter]] <- sample(1:n, cv.samp.size[iter], replace = FALSE)  
    } else {
      tmp <- c()
      for(j in 1:iter){
        tmp <- c(tmp, samp.cv[[j]])
      }
      samp.cv[[iter]] <- sample((1:n)[ - tmp], cv.samp.size[iter], replace = FALSE)
    } 
    Y.cv <- Y[ - samp.cv[[iter]]]
    Y.val <- Y[samp.cv[[iter]]]
    X.cv <- X[ - samp.cv[[iter]], ]
    X.val <- X[samp.cv[[iter]], ]
    out.cv[[iter]] <- list(X.cv = X.cv, Y.cv = Y.cv, X.val = X.val, Y.val = Y.val)
  }
  return(out.cv)
}