vechinv <- function(v) {
  # Inverse of vech function
  L <- vech2mat(v)
  L[upper.tri(L)] <- 0
  return(L)
}

A <- function(Q, r) {
  # Log-partition function of the multivariate normal distribution
  return(0.5*(t(r)%*%solve(Q)%*%r + log(det(2*pi*solve(Q)))))
}

tied_energy <- function(lambda, Z, Q_p, r_p, r_samples) {
  # Tied energy function for probit regression
  N <- nrow(Z)
  p <- ncol(Z)
  M <- nrow(r_samples)
  
  Q_L <- vechinv(lambda[1:(0.5*p*(p + 1))])
  diag(Q_L) <- softplus(diag(Q_L))
  Q <- Q_L%*%t(Q_L)
  r <- lambda[(0.5*p*(p + 1) + 1):(0.5*p*(p + 1) + p)]
  
  Q_dot <- Q_p + N*Q
  r_dot <- r_p + N*r
  
  Sigma_dot <- solve(Q_dot)
  mu_dot <- Sigma_dot%*%r_dot
  
  Sigma_dot_L <- t(chol(Sigma_dot))
  total <- 0
  
  for (n in 1:N) {
    subtotal <- 0
    
    for (m in 1:M) {
      q_sample <- Sigma_dot_L%*%r_samples[m, ] + mu_dot
      subtotal <- subtotal + pnorm(t(Z[n, ])%*%q_sample)/exp(-0.5*t(q_sample)%*%Q%*%q_sample + t(q_sample)%*%r)
    }
    
    total <- total + log(subtotal/M)
  }
  
  return(A(Q_p, r_p) - A(Q_dot, r_dot) - total)
}

bb_ep <- function(X, y, mu_beta, Sigma_beta,
                  M, epsilon, tau, r_init, Q_init,
                  min_iter, max_iter, thresh, verbose) {
  # Black-box EP for probit regression
  N <- nrow(X)
  p <- ncol(X)
  Z <- X*(2*y - 1)
  
  # Parameter initialisation
  Q <- Q_init
  r <- r_init
  
  Q_L <- t(chol(Q))
  diag(Q_L) <- softplusinv(diag(Q_L))
  lambda <- c(vech(Q_L), r)
  
  Q_p <- solve(Sigma_beta)
  r_p <- Q_p%*%mu_beta
  
  # Main EP loop
  for (iter in 1:max_iter) {
    if (verbose) print(paste0("---- Current iteration: ", iter, " ----"))
    
    r_samples <- matrix(rnorm(M*p), M)
    
    energy <- tied_energy(lambda, Z, Q_p, r_p, r_samples)
    energy_grad <- gradient(tied_energy, lambda, list(Z = Z, Q_p = Q_p, r_p = r_p, r_samples = r_samples))
    norm_energy_grad <- norm(energy_grad, "2")
    
    if (verbose) {
      print(paste0("Tied energy estimate:    ", energy))
      print(paste0("Norm of energy gradient: ", norm_energy_grad))
    }
    
    if (norm_energy_grad < thresh) {
      if (verbose) print("EP has converged; stopping EP")
      break
    } else {
      lambda <- lambda - min(epsilon, epsilon*tau/iter)*energy_grad
    }
  }
  
  # Returning in mean parameterisation
  Q_L <- vechinv(lambda[1:(0.5*p*(p + 1))])
  diag(Q_L) <- softplus(diag(Q_L))
  Q <- Q_L%*%t(Q_L)
  r <- lambda[(0.5*p*(p + 1) + 1):(0.5*p*(p + 1) + p)]
  
  Q_dot <- Q_p + N*Q
  r_dot <- r_p + N*r
  
  Sigma_dot <- solve(Q_dot)
  mu_dot <- Sigma_dot%*%r_dot
  
  return(list(mu = mu_dot, Sigma = Sigma_dot))
}
