// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(mvtnorm)]]
#include <mvtnormAPI.h>

// [[Rcpp::export]]
double pmvnorm_cppr(arma::vec& lower,
                    arma::vec& bound,
                    arma::vec& lowertrivec,
                    double abseps = 1e-3){
  
  int n = 2;
  int nu = 0;
  int maxpts = 25000;     // default in mvtnorm: 25000
  double releps = 0;      // default in mvtnorm: 0
  int rnd = 1;            // Get/PutRNGstate
  
  double* bound_ = bound.memptr();
  double* correlationMatrix = lowertrivec.memptr();
  double* lower_ = lower.memptr();
  int* infin = new int[n];
  double* delta = new double[n];
  
  for (int i = 0; i < n; ++i) {
    infin[i] = 2;
    if((lower_[i]==-1*arma::datum::inf) && (bound_[i]!=arma::datum::inf)) infin[i] = 0;
    if((lower_[i]!=-1*arma::datum::inf) && (bound_[i]==arma::datum::inf)) infin[i] = 1;
    if((lower_[i]==-1*arma::datum::inf) && (bound_[i]==arma::datum::inf)) infin[i] = -1;
    delta[i] = 0.0;
  }

//for (int i = 0; i < n; ++i) {
//  infin[i] = 2;
//  if(lower_[i] <-999) lower_[i] = -999.0;
//  if(bound_[i] > 999) bound_[i] = 999.0;
//  delta[i] = 0.0;
//}



  // return values
  double error;
  double value;
  int inform;
  
  mvtnorm_C_mvtdst(&n, &nu, lower_, bound_,
                   infin, correlationMatrix, delta,
                   &maxpts, &abseps, &releps,
                   &error, &value, &inform, &rnd);
                   // delete[] (lower);
                   delete[] (infin);
                   delete[] (delta);
                   
                   return value;
}