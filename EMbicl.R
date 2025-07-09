EMbicl<-function(theta,A,U,no.TH,pg,mX,Ind,G,P,Q,N,tol=10^-6,maxit=500,cl)
{
  # the function combines the E.step and the M.step until convergence
  # here the criterion of convergence is given by the difference
  # of likelihoods.
  #
  cat("|-------------|-------------|-------------|-------------|\n")
  cat("|     iter    |    imp      |     lik     |  linf-llik  |\n")
  cat("|-------------|-------------|-------------|-------------|\n")
  it <- 0
  likold <- -Inf
  likold <- -10^10
  dif <- Inf
  linf=rep(0,3)
  thetav <- parl2vr (theta,P,no.TH,G)
  A0<-A
  A<-A_update (theta,U,Ind,P,Q,A0,no.TH,G,cl)
  M <- M.step1(thetav,U,mX,Ind,P,A,no.TH,G,cl)
  
  thetav <- M$thetav
  pg <- M$pg
 
  
  theta <- parv2lr(list(thetav$param,thetav$paramn),A,P,G)
  
  Sigma=list()
  for (g in 1:G){
    Sigma[[g]]=A[[g]]%*%t(A[[g]])+theta$D[,,g];
  }
  theta$Sigma= Sigma
  theta$Mu=theta$Eta
  #  
  while((dif > tol) &  it < maxit) {
    it <- it+1
    U <- E.step1(theta,pg,Ind,G,cl)
    entr <- -sum(t(Ind[,6])%*%(U*log(U)))
    #print(system.time(M <- M.step1c(theta,U,mX,Ind,P,no.TH,G)))
    A0<-A
    A<-A_update (theta,U,Ind,P,Q,A0,no.TH,G,cl)
 print(A)
 
    M <- M.step1(thetav,U,mX,Ind,P,A,no.TH,G,cl)
    thetav <- M$thetav
    pg <- M$pg
    
    theta <- parv2lr(list(thetav$param,thetav$paramn),A,P,G)
    
    # theta$Eta=matrix(c( rep(0,P),
    #                     5, 5, 3, 3, 3, 1, 1, 1, 
    #                     
    #                     -5, -5, -5, -5, -3, -3, -3,-3),nrow=G,byrow=T)
    # 
    print(theta$Eta)
    Sigma=list()
    for (g in 1:G){
      Sigma[[g]]=A[[g]]%*%t(A[[g]])+theta$D[,,g];
    }
    theta$Sigma= Sigma
    theta$Mu=theta$Eta
   
   
    mon<-list(U=U,M=M,pg=pg,theta=theta)
    # save("mon",file="mon1.RData")
    #lik <- M$flik+entr+t(Ind[,6])%*%U%*%(log(pg))
    
    lik <- M$flik+entr+t(Ind[,6])%*%U%*%t(log(pg))
    dif<-(lik-likold)
    #dif <- (lik-likold)/abs(likold)
    #if (it==1){ dif <- (lik-likold)}else{dif <- (lik-likold)/abs(likold)}
    
    likold <- lik
    cat(sprintf("%11g",c(it,M$improv,lik,dif)),"\n", sep = " | ")
  }
  cat("|-------------|-------------|-------------|-------------|\n")
  CLC=-2*lik+2*sum(U*log(U))
  #G*P matrice A
  npar=length(thetav$param)+G-1+G*P
  ICL=-2*lik+2*sum(U*log(U))+npar*log(N)
  out <- list(lik=lik,U=U,pg=pg,Eta=theta$Eta,Sigma=theta$Sigma,D=theta$D,
              gamma=theta$ts,A=A,
              CLC=CLC,ICL=ICL)
  return(out)
}