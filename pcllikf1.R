pcllikf1 <- function(thetav,thetan,U,Ind,P,A,no.TH,G){
  library("Rcpp")
  sourceCpp("pmvnormCPPr.cpp")
  source("parv2lr.R")
  # 
  # pairwise complete log-likelihood
  # thetav vector parameters, pg prior probabilities, Ind cells 4 indices to 1
  # P no variables, no.TH no categories, G no groups, mX bivariate marginals 
  #
  # debug: appena entrati
 # cat("DEBUG: length(thetav)=", length(thetav), "\n")
  
  
  theta <- parv2lr(list(thetav,thetan),A,P,G)
  Sigma=list()
  for (g in 1:G){
    Sigma[[g]]=A[[g]]%*%t(A[[g]])+theta$D[,,g];
  }
  
  theta$Sigma= Sigma
  theta$Mu=theta$Eta
  
  # theta$Mu=matrix(c( rep(0,P),
  #                    5, 5, 3, 3, 3, 1, 1, 1, 
  #                    
  #                    -5, -5, -5, -5, -3, -3, -3,-3),nrow=G,byrow=T)
  rU <- dim(Ind)[1]
  pcll <- 0
  #l=matrix(0,rU,G)
  L <- matrix(0,2,2)
  for(g in 1:G){
    ds <- 1/sqrt(diag(theta$Sigma[[g]]))
    R <- theta$Sigma[[g]]*(ds%*%t(ds))
    for(i in 1:rU){
      L[1,1] <- theta$ts[[Ind[i,1]]][Ind[i,3]]
      L[1,2] <- theta$ts[[Ind[i,2]]][Ind[i,4]]
      L[2,1] <- theta$ts[[Ind[i,1]]][Ind[i,3]+1]
      L[2,2] <- theta$ts[[Ind[i,2]]][Ind[i,4]+1]
      Lg <- (L-matrix(theta$Mu[g,Ind[i,1:2]],2,2,byrow=TRUE))*matrix(ds[Ind[i,1:2]],2,2,byrow=TRUE)
      set.seed(i)
      prig <- pmvnorm_cppr(Lg[1,],Lg[2,],R[Ind[i,1],Ind[i,2]],1e-3)
      #prig <- pmvnorm(Lg[1,],Lg[2,],corr=R[Ind[i,1:2],Ind[i,1:2]],1e-3)
      #l[i,g]=prig
      #if(prig<.Machine$double.eps){prig <- .Machine$double.eps}
      prig[prig<0.000001]=0.000001
      #if(prig < 1e-7){prig <- 1e-7}
      #print(pcll)
      #print(U[i,g]*Ind[i,6]*log(prig))
      
      pc=U[i,g]*Ind[i,6]*log(prig)
      
      #pc[pc<log(.Machine$double.eps)]=log(.Machine$double.eps)
      pcll <- pcll+pc
      
      
    }
  }
  # debug: appena prima del return
 # cat("DEBUG: typeof(pcll)=", typeof(pcll), 
#      " length(pcll)=", length(pcll), "\n")
  return(as.numeric(sum(pcll)))
    
}