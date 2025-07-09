t3postProb1 <- function(X,theta,pg,all.pairs,Ind,G,cl){
  library("parallel")
  library("Rcpp")
  sourceCpp("pmvnormCPPr.cpp")
  source("t3int_post.R")
  # 
  # pairwise complete log-likelihood
  # thetav vector parameters, pg prior probabilities, Ind cells 4 indices to 1
  # P no variables, no.TH no categories, G no groups, mX bivariate marginals 
  #
  
  rU <- dim(Ind)[1]
  pcll <- 0
  l1=matrix(0,rU,G)
  Resp<- parLapply(cl,1:G,t3int_post,theta,pg,Ind)
  l1=matrix(unlist(Resp),rU,G)
  N=dim(X)[1]
  Uh=matrix(0,N,G)
  for(n in 1:N){
    resp=X[n,]
    L1=matrix(0,dim(all.pairs)[2],G)
    for(j in 1:dim(all.pairs)[2]){
      sel=c(all.pairs[1,j],all.pairs[2,j],resp[all.pairs[1,j]],resp[all.pairs[2,j]])
      idx=which(apply(Ind[,1:4], 1, function(x) return(all(x == sel))))
      #Uh[n,]=Uh[n,]*(l1[idx,])
      L1[j,]=log(l1[idx,])
    }
    #v=Uh[n,]
    #v[v<.Machine$double.eps]=.Machine$double.eps
    #Uh[n,]=v
    Uh[n,]=apply(L1,2,sum)
  }
  U=Uh
  ind=rep(0,N)
  for (n in 1:N){
    idL1=which(Uh[n,]<(-700))
    U[n,idL1]=-700
  }
  
  #U=exp(U)%*%diag(pg)
  # U=exp(-7.0e2*ind+Uh*(1-ind))%*%diag(pg)
  U=exp(U)*matrix(rep(pg,N),N,G,byrow = T)
  Uhat=t(apply(U,1,function(x) x/sum(x)))
  return(Uhat)
}
