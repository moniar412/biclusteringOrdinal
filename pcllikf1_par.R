pcllikf1_par <- function(thetav, thetan, U, Ind, P, A, no.TH, G, n_par) {
  source("pcllikf1.R")
  # se thetav è un singolo punto...
  if (length(thetav) == n_par) {
    return(pcllikf1(thetav, thetan, U, Ind, P, A, no.TH, G))
  }
  # se thetav è un batch di m punti
  if (length(thetav) %% n_par != 0) {
    stop("Lunghezza inattesa: ", length(thetav))
  }
  m <- length(thetav) / n_par
  mat <- matrix(thetav, nrow = m, byrow = TRUE)
  res <- numeric(m)
  for (i in seq_len(m)) {
    res[i] <- pcllikf1(mat[i, ], thetan, U, Ind, P, A, no.TH, G)
  }
  return(res)
}