log_int_0 <- function(a, b) {
  # Log of integral of pnorm(a + b*x)*dnorm(x) from -Inf to Inf
  # Source: A table of normal integrals (Owen, 1980)
  return(pnorm(a/sqrt(1 + b^2), log.p = T))
}

int_1_shift <- function(a, b, c) {
  # Integral of x*pnorm(a + b*x)*dnorm(x)*exp(c) from -Inf to Inf
  # Source: A table of normal integrals (Owen, 1980)
  return((b/sqrt(1 + b^2))*exp(dnorm(a/sqrt(1 + b^2), log = T) + c))
}

int_2_shift <- function(a, b, c) {
  # Integral of x^2*pnorm(a + b*x)*dnorm(x)*exp(c) from -Inf to Inf
  # Source: The explicit form of expectation propagation for a simple statistical model (Kim and Wand, 2016)
  return(exp(pnorm(a/sqrt(1 + b^2), log.p = T) + c) - (a*b^2/(1 + b^2)^(3/2))*exp(dnorm(a/sqrt(1 + b^2), log = T) + c))
}

ti <- function(mu, sigma_2) {
  # Tilted inference for probit regression
  sigma <- sqrt(sigma_2)
  
  log_tm_0 <- log_int_0(mu, sigma)
  tm_1 <- sigma*int_1_shift(mu, sigma, -log_tm_0) + mu
  tm_2 <- sigma_2*int_2_shift(mu, sigma, -log_tm_0) + 2*mu*tm_1 - mu^2
  
  return(list(mu = tm_1, sigma_2 = tm_2 - tm_1^2))
}

dep <- function(X, y, mu_beta, Sigma_beta,
                beta, r_init, Q_init,
                min_pass, max_pass, thresh, verbose, mem) {
  # Damped EP for probit regression
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
      
      ti_res <- ti(mu_c_star, Sigma_c_star)
      Q_h_star <- solve(ti_res$sigma_2)
      r_h_star <- Q_h_star%*%ti_res$mu
      
      Q_tilde <- (1 - beta)*Q_values[, , n] + beta*Z[n, ]%*%(Q_h_star - Q_c_star)%*%t(Z[n, ])
      r_tilde <- (1 - beta)*r_values[, n] + beta*Z[n, ]%*%(r_h_star - r_c_star)
      
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
