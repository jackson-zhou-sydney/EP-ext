int_0 <- function(a, b) {
  # Integral of pnorm(a + b*x)*dnorm(x) from -Inf to Inf
  # Source: A table of normal integrals (Owen, 1980)
  return(pnorm(a/sqrt(1 + b^2)))
}

int_1 <- function(a, b) {
  # Integral of x*pnorm(a + b*x)*dnorm(x) from -Inf to Inf
  # Source: A table of normal integrals (Owen, 1980)
  return((b/sqrt(1 + b^2))*dnorm(a/sqrt(1 + b^2)))
}

int_2 <- function(a, b) {
  # Integral of x^2*pnorm(a + b*x)*dnorm(x) from -Inf to Inf
  # Source: The explicit form of expectation propagation for a simple statistical model (Kim and Wand, 2016)
  return(pnorm(a/sqrt(1 + b^2)) - (a*b^2/(1 + b^2)^(3/2))*dnorm(a/sqrt(1 + b^2)))
}

ti <- function(mu, sigma_2) {
  # Tilted inference for probit regression
  sigma <- sqrt(sigma_2)
  
  tm_0 <- int_0(mu, sigma)
  tm_1 <- sigma*int_1(mu, sigma) + mu*tm_0
  tm_2 <- sigma_2*int_2(mu, sigma) + 2*mu*tm_1 - mu^2*tm_0
  
  mu_h <- tm_1/tm_0
  sigma_2_h <- tm_2/tm_0 - mu_h^2
  
  return(list(mu = mu_h, sigma_2 = sigma_2_h))
}

dl_ep <- function(X, y, mu_beta, Sigma_beta,
                  min_inner, max_inner, thresh_inner, outer_freq, r_init, Q_init,
                  min_pass, max_pass, thresh, verbose) {
  # Double-loop EP for probit regression
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
  
  Sigma_init <- solve(Q_init)
  mu_init <- Sigma_init%*%r_init
  
  # Parameter initialisation
  Q_values <- array(0, c(p, p, N))
  r_values <- matrix(0, p, N)
  
  Sigma_values <- array(0, c(p, p, N))
  mu_values <- matrix(0, p, N)
  
  Q_p <- solve(Sigma_beta)
  r_p <- Q_p%*%mu_beta
  
  Q_dot <- Q_p
  r_dot <- r_p
  
  for (n in 1:N) {
    Q_values[, , n] <- Q_init
    r_values[, n] <- r_init
    
    Sigma_values[, , n] <- Sigma_init
    mu_values[, n] <- mu_init
    
    Q_dot <- Q_dot + Q_values[, , n]
    r_dot <- r_dot + r_values[, n]
  }
  
  Sigma_dot <- solve(Q_dot)
  mu_dot <- Sigma_dot%*%r_dot
  
  # Main EP loop
  for (pass in 1:max_pass) {
    # Delta initialisation
    deltas_Q <- rep(0, N)
    deltas_r <- rep(0, N)
    
    if (verbose) print(paste0("---- Current pass: ", pass, " ----"))
    
    for (n in sample(1:N)) {
      Q_values_n_init <- Q_values[, , n]
      r_values_n_init <- r_values[, n]
      
      Q_c_outer <- Q_dot - Q_values[, , n]
      r_c_outer <- r_dot - r_values[, n]
      
      Q_dot_hat <- Q_dot
      r_dot_hat <- r_dot
      
      for (i in 1:max_inner) {
        Q_c <- Q_dot_hat - Q_values[, , n]
        r_c <- r_dot_hat - r_values[, n]
        
        Sigma_c <- solve(Q_c)
        mu_c <- Sigma_c%*%r_c
        
        Sigma_c_star <- t(Z[n, ])%*%Sigma_c%*%Z[n, ]
        mu_c_star <- t(Z[n, ])%*%mu_c
        
        Q_c_star <- solve(Sigma_c_star)
        r_c_star <- Q_c_star%*%mu_c_star
        
        ti_res <- ti(mu_c_star, Sigma_c_star)
        Q_h_star <- solve(ti_res$sigma_2)
        r_h_star <- Q_h_star%*%ti_res$mu
        
        Q_h <- Z[n, ]%*%(Q_h_star - Q_c_star)%*%t(Z[n, ]) + Q_c
        r_h <- Z[n, ]%*%(r_h_star - r_c_star) + r_c
        
        Sigma_h <- solve(Q_h)
        mu_h <- Sigma_h%*%r_h
        
        Sigma_tilde <- Sigma_values[, , n] + 1*(Sigma_h - Sigma_dot)
        mu_tilde <- mu_values[, n] + 1*(mu_h - mu_dot)
        
        Sigma_tilde_d <- Sigma_tilde - Sigma_values[, , n]
        mu_tilde_d <- mu_tilde - mu_values[, n]
        
        Q_tilde <- solve(Sigma_tilde)
        r_tilde <- Q_tilde%*%mu_tilde
        
        Q_tilde_d <- Q_tilde - Q_values[, , n]
        r_tilde_d <- r_tilde - r_values[, n]
        
        delta_Q_inner <- norm(Q_tilde_d, "F")
        delta_r_inner <- norm(r_tilde_d, "2")
        
        Q_dot <- Q_dot + Q_tilde_d
        r_dot <- r_dot + r_tilde_d
        
        Sigma_dot <- Sigma_dot + Sigma_tilde_d
        mu_dot <- mu_dot + mu_tilde_d
        
        Q_values[, , n] <- Q_tilde
        r_values[, n] <- r_tilde
        
        Sigma_values[, , n] <- Sigma_tilde
        mu_values[, n] <- mu_tilde
        
        if (i %% outer_freq == 0) {
          Q_dot_hat <- Q_c_outer + Q_values[, , n]
          r_dot_hat <- r_c_outer + r_values[, n]
        }
        
        if (i == 1) {
          bd_Q_inner <- delta_Q_inner
          bd_r_inner <- delta_r_inner
          next
        }
        
        if (!is.numeric(delta_Q_inner) || is.nan(delta_Q_inner) || !is.numeric(delta_r_inner) || is.nan(delta_r_inner)) browser()
        
        if (delta_Q_inner < thresh_inner*bd_Q_inner && delta_r_inner < thresh_inner*bd_r_inner && i >= min_inner) {
          break
        }
      }
      
      deltas_Q[n] <- norm(Q_values[, , n] - Q_values_n_init, "F")
      deltas_r[n] <- norm(r_values[, n] - r_values_n_init, "2")

      # Refresh mean parameters
      Sigma_dot <- solve(Q_dot)
      mu_dot <- Sigma_dot%*%r_dot
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
