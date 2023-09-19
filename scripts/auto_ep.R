ti_mc <- function(mu, sigma_2, M) {
  # Tilted inference for probit regression using Monte Carlo
  samples <- rnorm(M, mu, sqrt(sigma_2))
  
  tm_0 <- mean(sapply(samples, function(x) pnorm(x)))
  tm_1 <- mean(sapply(samples, function(x) x*pnorm(x)))
  tm_2 <- mean(sapply(samples, function(x) x^2*pnorm(x)))
  
  mu_h <- tm_1/tm_0
  sigma_2_h <- tm_2/tm_0 - mu_h^2
  
  return(list(mu = mu_h, sigma_2 = sigma_2_h))
}

gen_train <- function(X, y, mu_beta, Sigma_beta,
                      M, r_init, Q_init,
                      min_pass, max_pass, thresh, pat, verbose) {
  # Generate training data for automatic EP for probit regression
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
  wait <- 0
  
  train_mat <- matrix(NA, max_pass*N, 4)
  colnames(train_mat) <- c("mu", "sigma_2", "mu_h", "sigma_2_h")
  
  # Parameter initialisation
  Q_values <- array(0, c(p, p, N))
  r_values <- matrix(0, p, N)
  
  Q_p <- solve(Sigma_beta)
  r_p <- Q_p%*%mu_beta
  
  Q_dot <- Q_p
  r_dot <- r_p
  
  for (n in 1:N) {
    Q_values[, , n] <- Q_init
    r_values[, n] <- r_init
    
    Q_dot <- Q_dot + Q_values[, , n]
    r_dot <- r_dot + r_values[, n]
  }
  
  # Main EP loop
  for (pass in 1:max_pass) {
    # Delta initialisation
    deltas_Q <- rep(0, N)
    deltas_r <- rep(0, N)
    
    if (verbose) print(paste0("---- Current pass: ", pass, " ----"))
    
    for (n in sample(1:N)) {
      Q_c <- Q_dot - Q_values[, , n]
      r_c <- r_dot - r_values[, n]
      
      Sigma_c <- solve(Q_c)
      mu_c <- Sigma_c%*%r_c
      
      Sigma_c_star <- t(Z[n, ])%*%Sigma_c%*%Z[n, ]
      mu_c_star <- t(Z[n, ])%*%mu_c
      
      Q_c_star <- solve(Sigma_c_star)
      r_c_star <- Q_c_star%*%mu_c_star
      
      ti_mc_res <- ti_mc(mu_c_star, Sigma_c_star, M)
      train_mat[(pass - 1)*N + n, "mu"] <- mu_c_star
      train_mat[(pass - 1)*N + n, "sigma_2"] <- Sigma_c_star
      train_mat[(pass - 1)*N + n, "mu_h"] <- ti_mc_res$mu
      train_mat[(pass - 1)*N + n, "sigma_2_h"] <- ti_mc_res$sigma_2
      Q_h_star <- solve(ti_mc_res$sigma_2)
      r_h_star <- Q_h_star%*%ti_mc_res$mu
      
      Q_tilde <- Z[n, ]%*%(Q_h_star - Q_c_star)%*%t(Z[n, ])
      r_tilde <- Z[n, ]%*%(r_h_star - r_c_star)
      
      Q_tilde_d <- Q_tilde - Q_values[, , n]
      r_tilde_d <- r_tilde - r_values[, n]
      
      deltas_Q[n] <- norm(Q_tilde_d, "F")
      deltas_r[n] <- norm(r_tilde_d, "2")
      
      Q_dot <- Q_dot + Q_tilde_d
      r_dot <- r_dot + r_tilde_d
      
      Q_values[, , n] <- Q_tilde
      r_values[, n] <- r_tilde
    }
    
    if (pass == 1) {
      # Base maximum deltas
      bmd_Q <- max(deltas_Q)
      bmd_r <- max(deltas_r)
      
      # Minimum maximum deltas
      mmd_Q <- bmd_Q
      mmd_r <- bmd_r
      
      if (verbose) {
        print(paste0("Maximum delta for Q: ", bmd_Q))
        print(paste0("Maximum delta for r: ", bmd_r))
      }
      
      next
    }
    
    md_Q <- max(deltas_Q)
    md_r <- max(deltas_r)
    
    if (md_Q < mmd_Q || md_r < mmd_r) {
      mmd_Q <- min(md_Q, mmd_Q)
      mmd_r <- min(md_r, mmd_r)
      wait <- 0
    } else {
      wait <- wait + 1
    }
    
    if (verbose) {
      print(paste0("Maximum delta for Q: ", md_Q))
      print(paste0("Maximum delta for r: ", md_r))
      print(paste0("Waited for ", wait, " passes out of ", pat))
    }
    
    if (md_Q < thresh*bmd_Q && md_r < thresh*bmd_r && pass >= min_pass) {
      if (verbose) print("EP has converged; stopping EP")
      break
    } else if (wait == pat) {
      if (verbose) print("Exceeded patience; stopping EP")
      break
    }
  }
  
  # Returning training data
  return(na.omit(as.data.frame(train_mat)))
}

knn <- function(x, X_train, y_1_train, y_2_train, k) {
  # k-nearest neighbours with uncertainty
  x <- as.numeric(x)
  X_train <- as.matrix(X_train)
  
  dists <- sqrt(colSums((t(X_train) - x)^2))
  small_inds <- order(dists)[1:k]
  
  return(list(y_1 = mean(y_1_train[small_inds]),
              y_2 = mean(y_2_train[small_inds]),
              u = max(dists[small_inds])))
}

knn_u_max <- function(X, y_1, y_2, D_1_max, D_2_max, k, folds) {
  # Calculate maximum uncertainty using cross-validation
  N <- nrow(X)
  inds <- sample(rep(1:folds, ceiling(N/folds)))[1:N]
  error_mat <- matrix(0, N, 3)
  colnames(error_mat) <- c("D_1", "D_2", "u")
  n_e <- 1
  
  for (i in 1:folds) {
    X_train   <- X[inds != i, ]
    y_1_train <- y_1[inds != i]
    y_2_train <- y_2[inds != i]
    X_test   <- X[inds == i, ]
    y_1_test <- y_1[inds == i]
    y_2_test <- y_2[inds == i]
    
    for (j in 1:nrow(X_test)) {
      knn_res <- knn(X_test[j, ], X_train, y_1_train, y_2_train, k)
      error_mat[n_e, ] <- c(abs(knn_res$y_1 - y_1_test[j]),
                            abs(knn_res$y_2 - y_2_test[j]),
                            knn_res$u)
      n_e <- n_e + 1
    }
  }
  
  return(min(max(error_mat[, "u"][error_mat[, "D_1"] < D_1_max]),
             max(error_mat[, "u"][error_mat[, "D_2"] < D_2_max])))
}

auto_ep <- function(X, y, mu_beta, Sigma_beta,
                    k, train_df, u_max, M, r_init, Q_init,
                    min_pass, max_pass, thresh, pat, verbose, mem) {
  # Automatic EP for probit regression using k-nearest neighbours
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
  wait <- 0
  
  # Parameter initialisation
  Q_values <- array(0, c(p, p, N))
  r_values <- matrix(0, p, N)
  
  Q_p <- solve(Sigma_beta)
  r_p <- Q_p%*%mu_beta
  
  Q_dot <- Q_p
  r_dot <- r_p
  
  for (n in 1:N) {
    Q_values[, , n] <- Q_init
    r_values[, n] <- r_init
    
    Q_dot <- Q_dot + Q_values[, , n]
    r_dot <- r_dot + r_values[, n]
  }
  
  # Main EP loop
  for (pass in 1:max_pass) {
    # Delta initialisation
    deltas_Q <- rep(0, N)
    deltas_r <- rep(0, N)
    
    if (verbose) print(paste0("---- Current pass: ", pass, " ----"))
    
    for (n in sample(1:N)) {
      Q_c <- Q_dot - Q_values[, , n]
      r_c <- r_dot - r_values[, n]
      
      Sigma_c <- solve(Q_c)
      mu_c <- Sigma_c%*%r_c
      
      Sigma_c_star <- t(Z[n, ])%*%Sigma_c%*%Z[n, ]
      mu_c_star <- t(Z[n, ])%*%mu_c
      
      Q_c_star <- solve(Sigma_c_star)
      r_c_star <- Q_c_star%*%mu_c_star
      
      knn_res <- knn(c(mu_c_star, Sigma_c_star), train_df[, c("mu", "sigma_2")], train_df$mu_h, train_df$sigma_2_h, k)
      if (knn_res$u < u_max) {
        ti_res <- list(mu = knn_res$y_1, sigma_2 = knn_res$y_2)
      } else {
        ti_res <- ti_mc(mu_c_star, Sigma_c_star, M)
        train_df[nrow(train_df) + 1, ] <- c(mu_c_star, Sigma_c_star, ti_res$mu, ti_res$sigma_2)
      }
      Q_h_star <- solve(ti_res$sigma_2)
      r_h_star <- Q_h_star%*%ti_res$mu
      
      Q_star_tilde <- Q_h_star - Q_c_star
      r_star_tilde <- r_h_star - r_c_star
      
      if (Q_star_tilde < 0) {
        Q_tilde <- tryCatch(as.matrix(nearPD(Z[n, ]%*%Q_star_tilde%*%t(Z[n, ]))$mat), error = function(e) NA)
        r_tilde <- Z[n, ]%*%r_star_tilde
        if (!is.matrix(Q_tilde)) next
      } else {
        Q_tilde <- Z[n, ]%*%Q_star_tilde%*%t(Z[n, ])
        r_tilde <- Z[n, ]%*%r_star_tilde
      }
      
      Q_tilde_d <- Q_tilde - Q_values[, , n]
      r_tilde_d <- r_tilde - r_values[, n]
      
      deltas_Q[n] <- norm(Q_tilde_d, "F")
      deltas_r[n] <- norm(r_tilde_d, "2")
      
      Q_dot <- Q_dot + Q_tilde_d
      r_dot <- r_dot + r_tilde_d
      
      Q_values[, , n] <- Q_tilde
      r_values[, n] <- r_tilde
      
      if (mem) return(NA)
    }
    
    if (pass == 1) {
      # Base maximum deltas
      bmd_Q <- max(deltas_Q)
      bmd_r <- max(deltas_r)
      
      # Minimum maximum deltas
      mmd_Q <- bmd_Q
      mmd_r <- bmd_r
      
      if (verbose) {
        print(paste0("Maximum delta for Q: ", bmd_Q))
        print(paste0("Maximum delta for r: ", bmd_r))
      }
      
      next
    }
    
    md_Q <- max(deltas_Q)
    md_r <- max(deltas_r)
    
    if (md_Q < mmd_Q || md_r < mmd_r) {
      mmd_Q <- min(md_Q, mmd_Q)
      mmd_r <- min(md_r, mmd_r)
      wait <- 0
    } else {
      wait <- wait + 1
    }
    
    if (verbose) {
      print(paste0("Maximum delta for Q: ", md_Q))
      print(paste0("Maximum delta for r: ", md_r))
      print(paste0("Waited for ", wait, " passes out of ", pat))
    }
    
    if (md_Q < thresh*bmd_Q && md_r < thresh*bmd_r && pass >= min_pass) {
      if (verbose) print("EP has converged; stopping EP")
      break
    } else if (wait == pat) {
      if (verbose) print("Exceeded patience; stopping EP")
      break
    }
  }
  
  # Returning in mean parameterisation
  Sigma_dot <- solve(Q_dot)
  mu_dot <- Sigma_dot%*%r_dot
  
  return(list(mu = mu_dot, Sigma = Sigma_dot))
}
