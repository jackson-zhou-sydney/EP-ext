log_h <- function(x, mu, sigma_2) {
  # Log density of tilted distribution
  return(pnorm(x, log.p = T) + dnorm(x, mu, sqrt(sigma_2), log = T))
}

log_h_grad <- function(x, mu, sigma_2) {
  # Gradient of log density of tilted distribution
  return(sn::zeta(1, x) - (x - mu)/sigma_2)
}

ti_l_site <- function(mu, sigma_2) {
  # Approximate tilted inference for probit regression using Laplace's method (site)
  x_star <- optim(mu, log_h, log_h_grad, mu = mu, sigma_2 = sigma_2,
                  control = list(fnscale = -1), method = "BFGS")$par
  
  Q <- -sn::zeta(2, x_star)
  r <- Q*x_star
  
  return(list(r = r, Q = Q))
}

lp_site <- function(X, y, mu_beta, Sigma_beta,
                    r_init, Q_init,
                    min_pass, max_pass, thresh, verbose, mem) {
  # Laplace propagation using sites for probit regression
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
  
  # Parameter initialisation
  Q_values <- array(0, c(p, p, N + 1))
  r_values <- matrix(0, p, N + 1)
  
  Q_p <- solve(Sigma_beta)
  r_p <- Q_p%*%mu_beta
  
  Q_dot <- 0
  r_dot <- 0
  
  for (n in 1:(N + 1)) {
    Q_values[, , n] <- Q_init
    r_values[, n] <- r_init
    
    Q_dot <- Q_dot + Q_values[, , n]
    r_dot <- r_dot + r_values[, n]
  }
  
  # Main EP loop
  for (pass in 1:max_pass) {
    # Delta initialisation
    deltas_Q <- rep(0, N + 1)
    deltas_r <- rep(0, N + 1)
    
    if (verbose) print(paste0("---- Current pass: ", pass, " ----"))
    
    for (n in sample(1:(N + 1))) {
      Q_c <- Q_dot - Q_values[, , n]
      r_c <- r_dot - r_values[, n]
      
      if (n < N + 1) {
        Sigma_c <- solve(Q_c)
        mu_c <- Sigma_c%*%r_c
        
        Sigma_c_star <- t(Z[n, ])%*%Sigma_c%*%Z[n, ]
        mu_c_star <- t(Z[n, ])%*%mu_c
        
        Q_c_star <- solve(Sigma_c_star)
        r_c_star <- Q_c_star%*%mu_c_star
        
        ti_res <- ti_l_site(mu_c_star, Sigma_c_star)
        Q_star_tilde <- ti_res$Q
        r_star_tilde <- ti_res$r
        
        Q_tilde <- Z[n, ]%*%Q_star_tilde%*%t(Z[n, ])
        r_tilde <- Z[n, ]%*%r_star_tilde
      } else {
        mu_tilde <- solve(Q_c + Q_p)%*%(r_c + r_p)
        Q_tilde <- Q_p
        r_tilde <- Q_tilde%*%mu_tilde
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
  
  # Returning in mean parameterisation
  Sigma_dot <- solve(Q_dot)
  mu_dot <- Sigma_dot%*%r_dot
  
  return(list(mu = mu_dot, Sigma = Sigma_dot))
}
