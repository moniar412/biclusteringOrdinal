E.step1<-function(theta,pg,Ind,G,cl){
  #output: compute the E-step of the EM algorithm, i.e. we calculate the matrix 
  #with group membership probabilities for each subject. nrow=the total number of
  #pi^xi_ci (it is not equal to N!!!);ncol=the total number of groups, G.
  #input: theta is a list containing the parameters of each group (by row).
  #the parameter vector is composed of means,st.deviations,correlations and
  #the thresholds that are invariant across the groups (gamma).
  #th.idx is a vector keeping track to which variable the thresholds belongs to
  #X is the data
  #P is the number of variables 
  #R is the number of occasions
  #G is the number of groups
  #
  library("parallel")
  library("Rcpp")
  sourceCpp("pmvnormCPPr.cpp")
  source("int.R")
  rU <- dim(Ind)[1]
  Resp<- parLapply(cl,1:G,int,theta,pg,Ind)
  U=matrix(unlist(Resp),rU,G)
  U=U/apply(U,1,sum)
  return(U)
}