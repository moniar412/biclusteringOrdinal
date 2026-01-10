idx <- function(P,mX,no.TH){
  #
  # P no variables, mX bivariate marginals, no.TH no categories
  # output: var 1, var 2, cat 1, cat 2, table, cell freq
  #
  j <- 0
  Ind <- NULL
  for(p in 1:(P-1)){
    for(r in (p+1):P){
      j <- j+1
      for(c1 in 1:no.TH[p]){
        for(c2 in 1:no.TH[r]){
          fr <- mX[[j]][c1,c2]
          if(fr>0){Ind <- rbind(Ind,c(p,r,c1,c2,j,fr))}
        }
      }
    }
  }
  return(Ind)
} 