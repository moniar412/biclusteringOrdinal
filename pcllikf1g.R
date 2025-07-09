pcllikf1g <- function(theta,U,Ind,P,A,no.TH,G,g,cl){
  library("Rcpp")
  sourceCpp("pmvnormCPPr.cpp")
  source("parv2lr.R")
  # 
  # pairwise complete log-likelihood
  # thetav vector parameters, pg prior probabilities, Ind cells 4 indices to 1
  # P no variables, no.TH no categories, G no groups, mX bivariate marginals 
  #
  
  #theta <- parv2lr(list(thetav,thetan),P,G)
  pg <- (t(Ind[,6])%*%U)/sum(Ind[,6])
  #ord <- order(pg)
  #pg <- pg[ord]
  #theta$D<-theta$D[,,ord]
  #U=U[,ord]
  #Sigma=theta$Sigma[[g]]
  Sigma=A%*%t(A)+theta$D[,,g];
  #print(eigen(Sigma[[g]])$values)
  
  
  Mu=theta$Eta
  rU <- dim(Ind)[1]
  pcll <- 0
  #l=matrix(0,rU,G)
  L <- matrix(0,2,2)
  
  ds <- 1/sqrt(diag(Sigma))
  R <- Sigma*(ds%*%t(ds))
  for(i in 1:rU){
    L[1,1] <- theta$ts[[Ind[i,1]]][Ind[i,3]]
    L[1,2] <- theta$ts[[Ind[i,2]]][Ind[i,4]]
    L[2,1] <- theta$ts[[Ind[i,1]]][Ind[i,3]+1]
    L[2,2] <- theta$ts[[Ind[i,2]]][Ind[i,4]+1]
    Lg <- (L-matrix(Mu[g,Ind[i,1:2]],2,2,byrow=TRUE))*matrix(ds[Ind[i,1:2]],2,2,byrow=TRUE)
    set.seed(i)
    prig <- pmvnorm_cppr(Lg[1,],Lg[2,],R[Ind[i,1],Ind[i,2]],1e-3)
    prig[prig<0.000001]=0.000001
    #prig[prig<.Machine$double.eps]=.Machine$double.eps
    pc=U[i,g]*Ind[i,6]*log(prig)+Ind[i,6]*U[i,g]*log(pg[g])
    pcll <- pcll+pc
  }
  #entr <- -sum(t(Ind[,6])*(U[,g]*log(U[,g])))
  pcll<-pcll
  #+entr+sum(t(Ind[,6])*U[,g]*rep(log(pg[g]),length(U[,g])))
  #print(pcll)
  return(pcll)
}