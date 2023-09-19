log_h <- function(x, mu, sigma_2) {
  # Log density of tilted distribution
  return(pnorm(x, log.p = T) + dnorm(x, mu, sqrt(sigma_2), log = T))
}

log_h_grad <- function(x, mu, sigma_2) {
  # Gradient of log density of tilted distribution
  return(sn::zeta(1, x) - (x - mu)/sigma_2)
}

log_h_hess <- function(x, mu, sigma_2) {
  # Hessian of log density of tilted distribution
  return(sn::zeta(2, x) - 1/sigma_2)
}

ti_l_tilted <- function(mu, sigma_2) {
  # Approximate tilted inference for probit regression using Laplace's method (tilted)
  x_star <- optim(mu, log_h, log_h_grad, mu = mu, sigma_2 = sigma_2,
                  control = list(fnscale = -1), method = "BFGS")$par
  
  return(list(mu = x_star, sigma_2 = -solve(log_h_hess(x_star, mu, sigma_2))))
}

lp_tilted <- function(X, y, mu_beta, Sigma_beta,
                      r_init, Q_init,
                      min_pass, max_pass, thresh, verbose, mem) {
  # Laplace propagation using tilted distributions for probit regression
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
  
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
      
      ti_res <- ti_l_tilted(mu_c_star, Sigma_c_star)
      Q_h_star <- solve(ti_res$sigma_2)
      r_h_star <- Q_h_star%*%ti_res$mu
      
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
