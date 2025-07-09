t3int_post<-function(g,theta,pg,Ind){
  require("Rcpp")
  sourceCpp("pmvnormCPPr.cpp")
  ds <- 1/sqrt(diag(theta$Sigma[[g]]))
  R <- theta$Sigma[[g]]*(ds%*%t(ds))
  rU <- dim(Ind)[1]
  Ui=rep(0,rU)
  for(i in 1:rU){
    L <- matrix(0,2,2)
    L[1,1] <- theta$ts[[Ind[i,1]]][Ind[i,3]]
    L[1,2] <- theta$ts[[Ind[i,2]]][Ind[i,4]]
    L[2,1] <- theta$ts[[Ind[i,1]]][Ind[i,3]+1]
    L[2,2] <- theta$ts[[Ind[i,2]]][Ind[i,4]+1]
    Lg <- (L-matrix(theta$Mu[g,Ind[i,1:2]],2,2,byrow=TRUE))*matrix(ds[Ind[i,1:2]],2,2,byrow=TRUE)
    set.seed(i)
    uig <- pmvnorm_cppr(Lg[1,],Lg[2,],R[Ind[i,1],Ind[i,2]],1e-3)
    #if(uig<.Machine$double.eps){uig <- .Machine$double.eps}
    uig[uig< 1e-7]=1e-7
    #if(uig < 1e-7){uig <- 1e-7}
    Ui[i]<- uig
    
  }
  return(Ui)
}