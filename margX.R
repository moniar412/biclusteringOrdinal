margX <- function(X,no.TH){
  #
  # computes all the bivariate marginals
  #
  P <- dim(X)[2]
  mX <- list()
  i <- 0
  for(p in 1:(P-1)){
    x <- factor(X[,p],levels=1:no.TH[p])
    for(r in (p+1):P){
      i <- i+1
      y <- factor(X[,r],levels=1:no.TH[r])
      mX[[i]] <- table(x,y)
    }
  }
  return(mX)
}