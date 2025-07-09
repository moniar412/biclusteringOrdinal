# Estimate thresholds via ECDF (for each variable)
# Estimate thresholds via ECDF per variable
estimate_thresholds <- function(data) {
  apply(data, 2, function(x) {
    tab <- prop.table(table(x))
    cum_probs <- cumsum(tab)
    qnorm(cum_probs[-length(cum_probs)])  # remove 1 to avoid Inf
  }, simplify = FALSE)
}


# ------------------------------------
# Estimate latent means and variances
# ------------------------------------
compute_mean_var <- function(data, labels, thresholds_list, base_group = 1) {
  K <- length(unique(labels))
  P <- ncol(data)
  means <- matrix(0, nrow = K, ncol = P)
  vars   <- matrix(0, nrow = K, ncol = P)
  
  for (k in 1:K) {
    for (j in 1:P) {
      # Empirical probabilities for each level
      x <- data[labels == k, j]
      level_counts <- table(factor(x, levels = 1:(length(thresholds_list[[j]]) + 1)))
      probs <- prop.table(level_counts)
      
      # Breaks from estimated thresholds
      breaks <- c(-Inf, thresholds_list[[j]], Inf)
      
      # Midpoints for approximation
      mids <- sapply(1:(length(breaks)-1), function(i) {
        if (is.finite(breaks[i]) & is.finite(breaks[i+1])) {
          mean(c(breaks[i], breaks[i+1]))
        } else if (!is.finite(breaks[i])) {
          breaks[i+1] - 1
        } else {
          breaks[i] + 1
        }
      })
      
      means[k, j] <- sum(probs * mids)
      vars[k, j]   <- (sum(probs * (mids - means[k, j])^2))
    }
  }
  
  # Normalize base group to mean 0, sd 1
  for (j in 1:P) {
    means[, j] <- (means[, j] - means[base_group, j]) / vars[base_group, j]
    vars[, j]   <- vars[, j] / vars[base_group, j]
  }
  
  list(means = means, vars = vars)
}


compute_mean_var_ml <- function(data, labels, thresholds_list, base_group = 1) {
  K <- length(unique(labels))
  P <- ncol(data)
  mu    <- matrix(0, nrow = K, ncol = P)
  sigma2 <- matrix(0, nrow = K, ncol = P)
  
  for (k in 1:K) {
    idx_k <- which(labels == k)
    for (j in 1:P) {
      x_kj <- data[idx_k, j]
      # costruiamo i break dai thresholds
      breaks <- c(-Inf, thresholds_list[[j]], Inf)
      # conteggi osservati per categoria
      counts <- as.numeric(table(factor(x_kj, levels = 1:(length(breaks)-1))))
      # funzione di log‐verosimiglianza negativa
      nll <- function(par) {
        mu0   <- par[1]
        logsd <- par[2]
        sd0   <- exp(logsd)
        # probabilità teoriche di cadere in ogni categoria
        p_cat <- pnorm((breaks[-1] - mu0)/sd0) -
          pnorm((breaks[-length(breaks)] - mu0)/sd0)
        # evitiamo log(0)
        p_cat[p_cat < 1e-12] <- 1e-12
        -sum(counts * log(p_cat))
      }
      # stima via optim
      opt <- optim(c(0, 0), nll, method = "BFGS", hessian = FALSE)
      mu[k, j]    <- opt$par[1]
      sigma2[k, j] <- exp(2 * opt$par[2])
    }
  }
  
  # --- normalizzazione base_group su mu=0, var=1 ---
  for (j in 1:P) {
    mu[, j]    <- (mu[, j] - mu[base_group, j]) / sqrt(sigma2[base_group, j])
    sigma2[, j] <-  sigma2[, j]   /  sigma2[base_group, j]
  }
  
  list(means = mu, vars = sigma2)
}



library(polycor)

auto_collapse_ordinal_loop <- function(v1, v2, min_count = 5, max_iter = 10) {
  v1_new <- as.numeric(as.character(v1))
  v2_new <- as.numeric(as.character(v2))
  log_merge <- list(v1 = list(), v2 = list())
  
  for (i in 1:max_iter) {
    tbl <- table(v1_new, v2_new)
    mat <- as.matrix(tbl)
    
    bad_rows <- which(rowSums(mat) < min_count)
    bad_cols <- which(colSums(mat) < min_count)
    
    if (length(bad_rows) == 0 && length(bad_cols) == 0) break
    
    for (r in bad_rows) {
      current <- as.numeric(rownames(mat)[r])
      levels <- sort(unique(v1_new))
      pos <- which(levels == current)
      if (pos > 1) {
        new_val <- levels[pos - 1]
      } else if (pos < length(levels)) {
        new_val <- levels[pos + 1]
      } else {
        next
      }
      v1_new[v1_new == current] <- new_val
      log_merge$v1[[as.character(current)]] <- new_val
    }
    
    for (c in bad_cols) {
      current <- as.numeric(colnames(mat)[c])
      levels <- sort(unique(v2_new))
      pos <- which(levels == current)
      if (pos > 1) {
        new_val <- levels[pos - 1]
      } else if (pos < length(levels)) {
        new_val <- levels[pos + 1]
      } else {
        next
      }
      v2_new[v2_new == current] <- new_val
      log_merge$v2[[as.character(current)]] <- new_val
    }
  }
  
  return(list(
    v1_collapsed = factor(v1_new, ordered = TRUE),
    v2_collapsed = factor(v2_new, ordered = TRUE),
    log_merge = log_merge
  ))
}

auto_clean_dataset_polychoric <- function(data, min_count = 5, max_iter = 10) {
  var_names <- names(data)
  n <- ncol(data)
  data_clean <- data
  merge_log <- list()
  cor_matrix <- matrix(NA, nrow = n, ncol = n, dimnames = list(var_names, var_names))
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      v1 <- data_clean[[i]]
      v2 <- data_clean[[j]]
      
      result <- tryCatch({
        polycor::polychor(v1, v2)
      }, error = function(e) NA)
      
      if (!is.na(result)) {
        cor_matrix[i, j] <- cor_matrix[j, i] <- result
      } else {
        # Try collapsing and recompute
        collapsed <- auto_collapse_ordinal_loop(v1, v2, min_count, max_iter)
        cor_fixed <- tryCatch({
          polycor::polychor(collapsed$v1_collapsed, collapsed$v2_collapsed)
        }, error = function(e) NA)
        
        cor_matrix[i, j] <- cor_matrix[j, i] <- cor_fixed
        
        # Replace cleaned variables in the dataset
        data_clean[[i]] <- collapsed$v1_collapsed
        data_clean[[j]] <- collapsed$v2_collapsed
        
        # Log merges
        merge_log[[paste0(var_names[i], "_", var_names[j])]] <- collapsed$log_merge
      }
    }
  }
  
  diag(cor_matrix) <- 1
  return(list(data_clean = data_clean, cor_matrix = cor_matrix, merge_log = merge_log))
}

