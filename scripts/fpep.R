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

fpep <- function(X, y, mu_beta, Sigma_beta,
                 r_init, Q_init,
                 min_pass, max_pass, thresh, verbose) {
  # Fully-parallel EP for probit regression
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
  
  # Parameter initialisation
  Q_p <- solve(Sigma_beta)
  r_p <- Q_p%*%mu_beta
  
  foreach_res <- foreach(n = 1:N) %dopar% {
    Q_tilde <- Q_init
    r_tilde <- r_init
    
    Q_dot_local <- Q_tilde
    r_dot_local <- r_tilde
    
    return(list(Q_tilde = Q_tilde, 
                r_tilde = r_tilde,
                Q_dot_local = Q_dot_local, 
                r_dot_local = r_dot_local))
  }
  
  transpose_res <- transpose(foreach_res)
  
  Q_values <- abind(transpose_res$Q_tilde, along = 3)
  r_values <- do.call(cbind, transpose_res$r_tilde)
  
  Q_dot <- Q_p + Reduce("+", transpose_res$Q_dot_local)
  r_dot <- r_p + Reduce("+", transpose_res$r_dot_local)
  
  # Main EP loop
  for (pass in 1:max_pass) {
    if (verbose) print(paste0("---- Current pass: ", pass, " ----"))
    
    foreach_res <- foreach(n = 1:N, .export = c("int_0", "int_1", "int_2", "ti")) %dopar% {
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
      
      Q_tilde <- Z[n, ]%*%(Q_h_star - Q_c_star)%*%t(Z[n, ])
      r_tilde <- Z[n, ]%*%(r_h_star - r_c_star)
      
      Q_tilde_d <- Q_tilde - Q_values[, , n]
      r_tilde_d <- r_tilde - r_values[, n]
      
      deltas_Q_local <- norm(Q_tilde_d, "F")
      deltas_r_local <- norm(r_tilde_d, "2")
      
      Q_dot_local <- Q_tilde
      r_dot_local <- r_tilde
      
      return(list(Q_tilde = Q_tilde, 
                  r_tilde = r_tilde,
                  deltas_Q_local = deltas_Q_local,
                  deltas_r_local = deltas_r_local,
                  Q_dot_local = Q_dot_local, 
                  r_dot_local = r_dot_local))
    }
    
    transpose_res <- transpose(foreach_res)
    
    Q_values <- abind(transpose_res$Q_tilde, along = 3)
    r_values <- do.call(cbind, transpose_res$r_tilde)
    
    deltas_Q <- unlist(transpose_res$deltas_Q_local)
    deltas_r <- unlist(transpose_res$deltas_r_local)
    
    Q_dot <- Q_p + Reduce("+", transpose_res$Q_dot_local)
    r_dot <- r_p + Reduce("+", transpose_res$r_dot_local)
    
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
