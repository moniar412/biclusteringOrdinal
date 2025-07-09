A_update <- function(theta, U, Ind, P, Q, A, no.TH, G, cl) {
  
  # Funzione interna per spostare una riga da colonna sovrappopolata a colonna vuota
  move_one <- function(Ag, q_empty, theta, U, Ind, P, no.TH, G, g, cl) {
    over_cols <- which(colSums(Ag) > 2)
    if (length(over_cols) > 0) {
      scores <- sapply(over_cols, function(col) {
        row_1s <- which(Ag[, col] == 1)
        if (length(row_1s) > 0) {
          test_A <- Ag
          test_A[row_1s[1], col] <- 0
          pcllikf1g(theta, U, Ind, P, test_A, no.TH, G, g, cl)
        } else {
          -Inf
        }
      })
      best_col <- over_cols[which.min(scores)]
      row_to_move <- which(Ag[, best_col] == 1)[1]
      Ag[row_to_move, best_col] <- 0
      Ag[row_to_move, q_empty] <- 1
    } else {
      available_rows <- which(rowSums(Ag) == 0)
      if (length(available_rows) > 0) {
        Ag[sample(available_rows, 1), q_empty] <- 1
      }
    }
    return(Ag)
  }
  
  for (g in 1:G) {
    Ag <- matrix(0, nrow = P, ncol = Q[g])
    Ag[1, 1] <- 1  # vincolo: prima variabile assegnata al primo cluster
    
    for (p in 2:P) {
      scores <- rep(NA, Q[g])
      for (q in 1:Q[g]) {
        test_A <- Ag
        test_A[p, q] <- 1
        scores[q] <- pcllikf1g(theta, U, Ind, P, test_A, no.TH, G, g, cl)
      }
      Ag[p, which.max(scores)] <- 1
    }
    
    # Controlla colonne vuote e sposta
    empty_cols <- which(colSums(Ag) == 0)
    for (q in empty_cols) {
      Ag <- move_one(Ag, q, theta, U, Ind, P, no.TH, G, g, cl)
    }
    
    # Secondo controllo per colonne ancora con un solo 1 (non strettamente necessario)
    under_cols <- which(colSums(Ag) == 1)
    for (q in under_cols) {
      Ag <- move_one(Ag, q, theta, U, Ind, P, no.TH, G, g, cl)
    }
    
    if (any(is.na(Ag))) stop("Error: NA in A")
    
    A[[g]] <- Ag
  }
  A[[1]]=A[[1]]*1/sqrt(Q[1])

  return(A)
}
