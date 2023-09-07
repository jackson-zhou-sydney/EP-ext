ljl <- function(beta, X, y, mu_beta, Sigma_beta) {
  Z <- X*(2*y - 1)
  return(as.numeric(sum(pnorm(Z%*%beta, log.p = T)) - 0.5*t(beta - mu_beta)%*%solve(Sigma_beta)%*%(beta - mu_beta)))
}

ljl_grad <- function(beta, X, y, mu_beta, Sigma_beta) {
  Z <- X*(2*y - 1)
  return(as.vector(t(Z)%*%zeta(1, Z%*%beta) - solve(Sigma_beta)%*%(beta - mu_beta)))
}

ljl_hess <- function(beta, X, y, mu_beta, Sigma_beta) {
  Z <- X*(2*y - 1)
  return(t(Z)%*%(as.vector(zeta(2, Z%*%beta))*Z) - solve(Sigma_beta))
}

ti_mc <- function(mu, sigma_2, M) {
  # Tilted inference for probit regression using Monte Carlo
  samples <- rnorm(M, mu, sqrt(sigma_2))
  
  tm_0 <- mean(sapply(samples, function(x) pnorm(x)))
  tm_1 <- mean(sapply(samples, function(x) x*pnorm(x)))
  tm_2 <- mean(sapply(samples, function(x) x^2*pnorm(x)))
  
  h_mu <- tm_1/tm_0
  h_sigma_2 <- tm_2/tm_0 - h_mu^2
  
  return(list(mu = h_mu, sigma_2 = h_sigma_2))
}

gen_train <- function(X, y, mu_beta, Sigma_beta, maxit, N_train, M) {
  # Generate probit regression tilted inference training data
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
  
  train_mat <- matrix(0, N_train, 4)
  colnames(train_mat) <- c("mu", "sigma_2", "mu_h", "sigma_2_h")
  
  # laplace_mu <- optim(rep(0, p), ljl, ljl_grad, 
  #                     X = X, y = y, mu_beta = mu_beta, Sigma_beta = Sigma_beta,
  #                     method = "BFGS", control = list(fnscale = -1, maxit = maxit))$par
  # laplace_Sigma <- -solve(ljl_hess(laplace_mu, X, y, mu_beta, Sigma_beta))
  
  for (n_t in 1:N_train) {
    n <- sample(1:N, 1)
    
    # sigma_2 <- t(Z[n, ])%*%rinvwishart(1, 100 + p, laplace_Sigma/10, checkSymmetry = F)[, , 1]%*%Z[n, ]
    # mu <- t(Z[n, ])%*%as.vector(rmvnorm(1, laplace_mu, 0.1*laplace_Sigma))
    sigma_2 <- rexp(1, rate = 2)
    mu <- rnorm(1, sd = 2)
    ti_mc_res <- ti_mc(mu, sigma_2, M)
    
    train_mat[n_t, "mu"] <- mu
    train_mat[n_t, "sigma_2"] <- sigma_2
    train_mat[n_t, "mu_h"] <- ti_mc_res$mu
    train_mat[n_t, "sigma_2_h"] <- ti_mc_res$sigma_2
  }
  
  return(na.omit(as.data.frame(train_mat)))
}

gen_train <- function(X, y, mu_beta, Sigma_beta,
                      M, r_init, Q_init,
                      min_pass, max_pass, thresh, verbose) {
  # Standard EP for probit regression
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
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
      
      # ti_mc_res <- ti_mc(mu_c_star, Sigma_c_star, M)
      ti_mc_res <- ti(mu_c_star, Sigma_c_star)
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
      deltas_r[n] <- norm(r_tilde_d, "F")
      
      Q_dot <- Q_dot + Q_tilde_d
      r_dot <- r_dot + r_tilde_d
      
      Q_values[, , n] <- Q_tilde
      r_values[, n] <- r_tilde
    }
    
    if (pass == 1) {
      # Base maximum deltas
      bmd_Q <- max(deltas_Q)
      bmd_r <- max(deltas_r)
      
      if (verbose) {
        print(paste0("Maximum delta for Q: ", bmd_Q))
        print(paste0("Maximum delta for r: ", bmd_r))
      }
      
      next
    }
    
    md_Q <- max(deltas_Q)
    md_r <- max(deltas_r)
    
    if (verbose) {
      print(paste0("Maximum delta for Q: ", md_Q))
      print(paste0("Maximum delta for r: ", md_r))
    }
    
    if (md_Q < thresh*bmd_Q && md_r < thresh*bmd_r && pass >= min_pass) {
      if (verbose) print("EP has converged; stopping EP")
      break
    }
  }
  
  # Returning training data
  
  return(na.omit(as.data.frame(train_mat)))
}

auto_ep <- function(X, y, mu_beta, Sigma_beta,
                    train_df, r_init, Q_init,
                    min_pass, max_pass, thresh, pat, verbose) {
  # Automatic EP for probit regression
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
  wait <- 0
  
  # Fitting machine learning model
  if (verbose) print("---- Fitting model ----")
  ti_rf_mu_h <- randomForest(mu_h ~ mu + sigma_2, data = train_df)
  ti_rf_sigma_2_h <- randomForest(sigma_2_h ~ mu + sigma_2, data = train_df)
  
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
      
      new_data <- data.frame(mu = mu_c_star, sigma_2 = Sigma_c_star)
      ti_res <- list(mu = predict(ti_rf_mu_h, new_data)[[1]], 
                     sigma_2 = predict(ti_rf_sigma_2_h, new_data)[[1]])
      ti_res_old <- ti(mu_c_star, Sigma_c_star)
      if (abs(ti_res$mu - ti_res_old$mu) > 0) ti_res <- ti_res_old
      Q_h_star <- solve(ti_res$sigma_2)
      r_h_star <- Q_h_star%*%ti_res$mu
      
      Q_tilde <- Z[n, ]%*%(Q_h_star - Q_c_star)%*%t(Z[n, ])
      r_tilde <- Z[n, ]%*%(r_h_star - r_c_star)
      
      Q_tilde_d <- Q_tilde - Q_values[, , n]
      r_tilde_d <- r_tilde - r_values[, n]
      
      deltas_Q[n] <- norm(Q_tilde_d, "F")
      deltas_r[n] <- norm(r_tilde_d, "F")
      
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
  
  # Returning in mean parameterisation
  Sigma_dot <- solve(Q_dot)
  mu_dot <- Sigma_dot%*%r_dot
  
  return(list(mu = mu_dot, Sigma = Sigma_dot))
}
