parv2lr <- function(param,A,P,G){
  #
  # from parameter vector to parameter list
  #
  mon1 <- param[[1]]
  mon2 <- param[[2]]
  # means
  eta <- rbind(rep(0,P),matrix(mon1[mon2=="Eta"],G-1,P))
  # variances
  d<-matrix(mon1[mon2=="d"],G-1,P,byrow=TRUE)
  D=array(0,c(P,P,G))
  #diag(D[,,1])<-1
  for (g in 1:(G-1)){
    #D[1,1,g+1]<-1
    diag(D[1:P,1:P,g+1])<-exp(d[g,])
  }
  D[,,1]=diag(1-diag(A[[1]]%*%t(A[[1]])))
  
  ts <- list()
  for(p in 1:P){
    a <- mon1[mon2==paste("ts",sep="",p)]
    ts[[p]] <- c(-Inf,c(a[1],a[1]+cumsum(exp(a[-1]))),Inf)
  }
  return(list(Eta=eta,D=D,ts=ts))
}