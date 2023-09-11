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

sep <- function(X, y, mu_beta, Sigma_beta,
                r_init, Q_init,
                min_pass, max_pass, thresh, verbose) {
  # Stochastic EP for probit regression
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
  
  # Parameter initialisation
  Q <- Q_init
  r <- r_init
  
  Q_p <- solve(Sigma_beta)
  r_p <- Q_p%*%mu_beta
  
  Q_dot <- Q_p + N*Q
  r_dot <- r_p + N*r
  
  # Main EP loop
  for (pass in 1:max_pass) {
    # Delta initialisation
    deltas_Q <- rep(0, N)
    deltas_r <- rep(0, N)
    
    if (verbose) print(paste0("---- Current pass: ", pass, " ----"))
    
    for (n in sample(1:N)) {
      Q_c <- Q_dot - Q
      r_c <- r_dot - r
      
      Sigma_c <- solve(Q_c)
      mu_c <- Sigma_c%*%r_c
      
      Sigma_c_star <- t(Z[n, ])%*%Sigma_c%*%Z[n, ]
      mu_c_star <- t(Z[n, ])%*%mu_c
      
      Q_c_star <- solve(Sigma_c_star)
      r_c_star <- Q_c_star%*%mu_c_star
      
      ti_res <- ti(mu_c_star, Sigma_c_star)
      Q_h_star <- solve(ti_res$sigma_2)
      r_h_star <- Q_h_star%*%ti_res$mu
      
      Q_tilde <- (1 - 1/N)*Q + (1/N)*Z[n, ]%*%(Q_h_star - Q_c_star)%*%t(Z[n, ])
      r_tilde <- (1 - 1/N)*r + (1/N)*Z[n, ]%*%(r_h_star - r_c_star)
      
      Q_tilde_d <- Q_tilde - Q
      r_tilde_d <- r_tilde - r
      
      deltas_Q[n] <- norm(Q_tilde_d, "F")
      deltas_r[n] <- norm(r_tilde_d, "2")
      
      Q_dot <- Q_dot + N*Q_tilde_d
      r_dot <- r_dot + N*r_tilde_d
      
      Q <- Q_tilde
      r <- r_tilde
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
