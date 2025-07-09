M.step1 <- function(thetav0,U,mX,Ind,P,A,no.TH,G,cl){
  #
  # M step
  #
  #
  # update pg
  pg <- (t(Ind[,6])%*%U)/sum(Ind[,6])
  #
  # update theta
  npar <- length(thetav0$param)
  oldvalue <- pcllikf1(thetav0$param,thetav0$paramn,U,Ind,P,A,no.TH,G)
  # set the number of processor cores
  setDefaultCluster(cl=cl) # set 'cl' as default cluster
  #clusterExport(cl,c("parv2l", "pmvnorm","pmvnorm_cpp","pnorm","triangl"))
  #
  
  ctrl <- 0
  while(ctrl==0){
    ctrl <- 1
    op <- optimParallel(par=thetav0$param,fn=pcllikf1,method="L-BFGS-B",
                        lower=-Inf,upper=Inf,thetan=thetav0$paramn,U=U,Ind=Ind,
                        P=P,A=A,no.TH=no.TH,G=G,
                        control=list(fnscale=-1))
                        #control=list(fnscale=-1,maxit=round(npar/4,0)))
    
    
    
    if((op$value-oldvalue)<0){ctrl <- 0; thetav0$param <- op$par}
  }
  
  thetav0$param <- op$par
  out <- list(thetav=thetav0,pg=pg,flik=op$value,improv=op$value-oldvalue)
  return(out)
}