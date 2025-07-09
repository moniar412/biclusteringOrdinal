simH1=function(N,P,Q,G,no.TH,nk,st1,cl){
  #medie e covarianze a blocchi
  
  all.pairs=combn(1:(P),2)
  
  
 
  
  #identification
  
  eta=matrix(c( rep(0,P),
                5, 5, 3, 3, 3, -1, -1, -1, 
                
                -3, -3,-3, -5, -5, -5,5,5),nrow=G,byrow=T)
  
    
 # eta=eta-matrix(rep(eta[1,],G),ncol=P,byrow=T)
#  eta[2,]=rep(4.5,8)
 # eta[3,]=rep(-4,8)
  A=list()
  A[[1]]=matrix(c(   1, 0,
                     1, 0,
                     1, 0,
                     1, 0,
                     0, 1,
                     0, 1,
                     0, 1,
                     0, 1),ncol=Q[1],byrow=T) 
  
  
  A[[1]]=A[[1]]*1/sqrt(Q[1])
  A[[2]]=matrix(c( 1, 0,0,
                   1, 0,0,
                   0, 1,0,
                   0, 1,0,
                   0, 1,0,
                   0, 0,1,
                   0, 0,1,
                   0, 0,1
  ),ncol=Q[2],byrow=T)
   A[[3]]=matrix(c(  1, 0, 0,
                     1, 0, 0,
                     1, 0, 0,
                     0, 1, 0,
                     0, 1, 0,
                     0, 1, 0,
                     0, 0, 1,
                     0, 0, 1),ncol=Q[3],byrow=T)
   
    
  
  D=array(0,c(P,P,G))
  D[,,1]=diag(1-diag(A[[1]]%*%t(A[[1]])))
  D[,,2]=diag(P)
  D[,,3]=diag(P)
  
  
  #D[5:8,5:8,1]=diag(4)*1.3
 # D[,,2]=diag(P)*seq(0.8,2,(2-0.8)/(P-1))
  D[1:2,1:2,2]=diag(2)*1
  D[3:5,3:5,2]=diag(3)*1.5
  D[6:8,6:8,2]=diag(3)*2
  
  
  #D[,,3]=diag(P)*seq(1.5,2.5,(2.5-1.5)/(P-1))
  D[1:3,1:3,3]=diag(3)*1
  D[4:6,4:6,3]=diag(3)*2
  D[7:8,7:8,3]=diag(2)*2.5
 
  
  
  
  
  
  # D[1,1,1]<-1
  # D[1,1,2]<-1
  # D[1,1,3]<-1
  
  Sigma=list()
  
  for (g in 1:G){
    Sigma[[g]]=A[[g]]%*%t(A[[g]])+D[,,g];
  }
  pg=c(.20, .30, .5)
  
  no.TH=rep(5,P)
  #cluster R, cluster C
  ARI=array(NA,dim=c(nk,G+1,length(N)))
  #Mean,Pg,D,ts
  Eu=array(NA,dim=c(nk,4,length(N)))
  for(n in 1:length(N)){
    for(rep in 1:nk){
      st=st1*nk+1+rep;
      set.seed(st)
      
      outx=randmixBicl(N[n],eta,Sigma,pg,no.TH,P,G,st)
      mX <- margX(outx$X,no.TH)
      Ind <- idx(P,mX,no.TH)  
      
      
      st0=stPoint(outx$X,outx$Y,G,Q,no.TH,st)
      theta=list()
      theta$Eta=st0$eta0
      theta$A=st0$A0
      theta$D=st0$D0
      theta$pg=st0$pg0
      theta$Sigma=st0$Sigma0
      theta$ts=st0$gamma0
      theta$Mu=theta$Eta
      theta$D=D
      Sigma0=list()
      for (g in 1:G){
        Sigma0[[g]]=A[[g]]%*%t(A[[g]])+theta$D[,,g];
      }
      theta$Sigma= Sigma0
      
      U=E.step1(theta,theta$pg,Ind,G,cl)
      #out=EMbicl(theta,theta$A,U,no.TH,theta$pg,mX,Ind,G,P,Q,N[n],tol=10^-1,maxit=500,cl)
      
     
      out=EMbicl(theta,theta$A,U,no.TH,theta$pg,mX,Ind,G,P,Q,N[n],tol=10^-3,maxit=500,cl)
      
      thetah=list()
      thetah$Eta=out$Eta
      thetah$A=out$A
      thetah$D=out$D
      thetah$pg=out$pg
      thetah$Sigma=out$Sigma
      thetah$ts=out$gamma
      thetah$Mu=out$Eta
      #classification
      Uhat=t3postProb1(outx$X,thetah,out$pg,all.pairs,Ind,G,cl)
      Hcl=apply(Uhat,1,which.max)
      #ARI
      ARI[rep,1,n]=adjustedRandIndex(outx$tcl,Hcl)
      for (g in 1:G){
        ARI[rep,g+1,n]=adjustedRandIndex(A[[g]],round(thetah$A[[g]]))
      }
      #computation of Euclidean squared distance between the estimates
      #and the true parameter values
      Eu[rep,1,n]=mean((thetah$Eta[2:G,]-eta[2:G,])^2)
      Eu[rep,2,n]=mean((thetah$pg[1:(G-1)]-pg[1:(G-1)])^2)
      Eu[rep,3,n]= mean((apply(thetah$D[,,2:G],3,diag)-apply(D[,,2:G],3,diag))^2)
      Eu[rep,4,n]=mean((unlist(thetah$ts)-unlist(outx$ts))^2,na.rm=T)
      
      write(ARI[rep,,n],file=paste(N[n],"N", "ARI_1.txt",sep=""),ncolumns=4,append=TRUE)
      write(Eu[rep,,n], file=paste(N[n],"N","ED_1.txt",sep=""),ncolumns=4,append=TRUE)
    }
  }
  out=list(ARI=ARI,Eu=Eu)
  return(out)
}