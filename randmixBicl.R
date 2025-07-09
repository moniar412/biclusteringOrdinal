randmixBicl<-function(N,Mu,Sigma,pg,no.TH,P,G,st)
{
  # output: random generation of a multivariate Gaussian mixture
  # input: N is the sample size
  # mu is a list of length G of mean vectors
  # Sigma is a list of length G of variance/covariance matrices
  # pg is a vector containing the probabilities associated with each group
  set.seed(st)
  n.samp <- sample(1:G,N,replace=TRUE,prob=pg)
  n.pg <- as.numeric(table(n.samp))
  cn.pg <- c(0,cumsum(n.pg))
  X <- matrix(0,N,P)
  for(g in 1:G){
    id <- cn.pg[g]+(1:n.pg[g])
    X[id,] <- rmvnorm(n=n.pg[g],mean=Mu[g,],sigma=Sigma[[g]])
  }
  ts <- list()
  for(p in 1:P){
    #st <- (max(X[,p])-min(X[,p]))/no.TH[p]
    #b <- seq(min(X[,p]),max(X[,p]),st)
    b<-quantile(X[,p], 1/(no.TH[p]+1)*1:(no.TH[p]+1))
    b[1] <- -Inf; b[no.TH[p]+1] <- Inf
    ts[[p]] <- b  
  }
  Y=X
  
  for(p in 1:(P)){
    X[,p] <- as.numeric(cut(X[,p],breaks=ts[[p]]))
  }
  tcl <- rep(1:G,times=n.pg)
  return(list(X=X,tcl=tcl,ts=ts,Y=Y))      
}