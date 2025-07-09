parl2vr <- function(theta,P,no.TH,G){
  #
  # from parameter list to vector: 3-way case
  #
  param <- c(theta$Eta[2:G,])
  paramn <- rep("Eta",(G-1)*P)
  #d=diag(D[,,1])
  d=0
  for (g in 2:G){
    d=c(d,diag(theta$D[1:P,1:P,g]))
    
  }
  d=d[-1]
  param <- c(param,log(d))
  paramn <- c(paramn,rep("d",(P)*(G-1)))
 
   for(p in 1:P){
    ts <- theta$ts[[p]][-c(1,no.TH[p]+1)]
    param <- c(param,c(ts[1],log(diff(ts))))
    paramn <- c(paramn,rep(paste("ts",sep="",p),no.TH[p]-1))
  }
  return(list(param=param,paramn=paramn))
}

