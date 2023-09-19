
# Numerical experiments

source("setup.R")
load("data/o_rings.RData")
load("results/mcmc.RData")
method <- "auto_ep"
set.seed(1)

init_list <- list()

for (i in 1:n_init) {
  init_list[[i]] <- list(r_init = rnorm(p), Q_init = diag(p))
}

l1_df <- data.frame(method = character(),
                    init = double(),
                    j = double(),
                    l1 = double())

other_df <- data.frame(method = character(),
                       init = double(),
                       time = double(),
                       mem = double())

if (method == "ep") {
  # Standard EP
  source("scripts/ep.R")
  
  for (i in 1:n_init) {
    print(paste0("Testing initialisation ", i))
    
    r_init <- init_list[[i]]$r_init
    Q_init <- init_list[[i]]$Q_init
    res <- tryCatch({set.seed(1); ep(X, y, mu_beta, Sigma_beta, 
                                     r_init, Q_init, 
                                     min_pass, max_pass, thresh, verbose, F)},
                    error = function(e) NA)
    if (is.list(res)) {
      for (j in 1:p) {
        values <- dnorm(grid_points[j, ], res$mu[j], sqrt(res$Sigma[j, j]))
        l1 <- 1 - 0.5*trapz(grid_points[j, ], abs(mcmc_values[j, ] - values))
        l1_df <- l1_df |> add_row(method = "ep", init = i, j = j, l1 = l1)
      }
      
      mark_res_1 <- mark({set.seed(1); ep(X, y, mu_beta, Sigma_beta, 
                                          r_init, Q_init, 
                                          min_pass, max_pass, thresh, verbose, F)})
      
      mark_res_2 <- mark({set.seed(1); ep(X, y, mu_beta, Sigma_beta, 
                                          r_init, Q_init, 
                                          min_pass, 1, thresh, verbose, T)})
      
      time <- as.numeric(mark_res_1$median)
      mem <- as.numeric(mark_res_2$mem_alloc)
      
      other_df <- other_df |> add_row(method = "ep", init = i, time = time, mem = mem)
    } else {
      for (j in 1:p) {
        l1_df <- l1_df |> add_row(method = "ep", init = i, j = j, l1 = NA)
      }
      
      other_df <- other_df |> add_row(method = "ep", init = i, time = NA, mem = NA)
    }
  }
  
  save(l1_df, other_df, file = "results/ep.RData")
} else if (method == "fpep") {
  # Distributed computation
  source("scripts/fpep.R")
  registerDoParallel(cores = fpep_cores)
  
  for (i in 1:n_init) {
    print(paste0("Testing initialisation ", i))
    
    r_init <- init_list[[i]]$r_init
    Q_init <- init_list[[i]]$Q_init
    res <- tryCatch({set.seed(1); fpep(X, y, mu_beta, Sigma_beta, 
                                       r_init, Q_init, 
                                       min_pass, max_pass, thresh, verbose, F)},
                    error = function(e) NA)
    if (is.list(res)) {
      for (j in 1:p) {
        values <- dnorm(grid_points[j, ], res$mu[j], sqrt(res$Sigma[j, j]))
        l1 <- 1 - 0.5*trapz(grid_points[j, ], abs(mcmc_values[j, ] - values))
        l1_df <- l1_df |> add_row(method = "fpep", init = i, j = j, l1 = l1)
      }
      
      mark_res_1 <- mark({set.seed(1); fpep(X, y, mu_beta, Sigma_beta, 
                                            r_init, Q_init, 
                                            min_pass, max_pass, thresh, verbose, F)})
      
      mark_res_2 <- mark({set.seed(1); fpep(X, y, mu_beta, Sigma_beta, 
                                            r_init, Q_init, 
                                            min_pass, 1, thresh, verbose, T)})
      
      time <- as.numeric(mark_res_1$median)
      mem <- as.numeric(mark_res_2$mem_alloc)/(N/fpep_cores)
      
      other_df <- other_df |> add_row(method = "fpep", init = i, time = time, mem = mem)
    } else {
      for (j in 1:p) {
        l1_df <- l1_df |> add_row(method = "fpep", init = i, j = j, l1 = NA)
      }
      
      other_df <- other_df |> add_row(method = "fpep", init = i, time = NA, mem = NA)
    }
  }
  
  save(l1_df, other_df, file = "results/fpep.RData")
  registerDoSEQ()
} else if (method == "lp_site") {
  # Laplace propagation using sites
  source("scripts/lp_site.R")
  
  for (i in 1:n_init) {
    print(paste0("Testing initialisation ", i))
    
    r_init <- init_list[[i]]$r_init
    Q_init <- init_list[[i]]$Q_init
    res <- tryCatch({set.seed(1); lp_site(X, y, mu_beta, Sigma_beta, 
                                          r_init, Q_init, 
                                          min_pass, max_pass, thresh, verbose, F)},
                    error = function(e) NA)
    if (is.list(res)) {
      for (j in 1:p) {
        values <- dnorm(grid_points[j, ], res$mu[j], sqrt(res$Sigma[j, j]))
        l1 <- 1 - 0.5*trapz(grid_points[j, ], abs(mcmc_values[j, ] - values))
        l1_df <- l1_df |> add_row(method = "lp_site", init = i, j = j, l1 = l1)
      }
      
      mark_res_1 <- mark({set.seed(1); lp_site(X, y, mu_beta, Sigma_beta, 
                                               r_init, Q_init, 
                                               min_pass, max_pass, thresh, verbose, F)})
      
      mark_res_2 <- mark({set.seed(1); lp_site(X, y, mu_beta, Sigma_beta, 
                                               r_init, Q_init, 
                                               min_pass, 1, thresh, verbose, T)})
      
      time <- as.numeric(mark_res_1$median)
      mem <- as.numeric(mark_res_2$mem_alloc)
      
      other_df <- other_df |> add_row(method = "lp_site", init = i, time = time, mem = mem)
    } else {
      for (j in 1:p) {
        l1_df <- l1_df |> add_row(method = "lp_site", init = i, j = j, l1 = NA)
      }
      
      other_df <- other_df |> add_row(method = "lp_site", init = i, time = NA, mem = NA)
    }
  }
  
  save(l1_df, other_df, file = "results/lp_site.RData")
} else if (method == "lp_tilted") {
  # Laplace propagation using tilted distributions
  source("scripts/lp_tilted.R")
  
  for (i in 1:n_init) {
    print(paste0("Testing initialisation ", i))
    
    r_init <- init_list[[i]]$r_init
    Q_init <- init_list[[i]]$Q_init
    res <- tryCatch({set.seed(1); lp_tilted(X, y, mu_beta, Sigma_beta, 
                                            r_init, Q_init, 
                                            min_pass, max_pass, thresh, verbose, F)},
                    error = function(e) NA)
    if (is.list(res)) {
      for (j in 1:p) {
        values <- dnorm(grid_points[j, ], res$mu[j], sqrt(res$Sigma[j, j]))
        l1 <- 1 - 0.5*trapz(grid_points[j, ], abs(mcmc_values[j, ] - values))
        l1_df <- l1_df |> add_row(method = "lp_tilted", init = i, j = j, l1 = l1)
      }
      
      mark_res_1 <- mark({set.seed(1); lp_tilted(X, y, mu_beta, Sigma_beta, 
                                                 r_init, Q_init, 
                                                 min_pass, max_pass, thresh, verbose, F)})
      
      mark_res_2 <- mark({set.seed(1); lp_tilted(X, y, mu_beta, Sigma_beta, 
                                                 r_init, Q_init, 
                                                 min_pass, 1, thresh, verbose, T)})
      
      time <- as.numeric(mark_res_1$median)
      mem <- as.numeric(mark_res_2$mem_alloc)
      
      other_df <- other_df |> add_row(method = "lp_tilted", init = i, time = time, mem = mem)
    } else {
      for (j in 1:p) {
        l1_df <- l1_df |> add_row(method = "lp_tilted", init = i, j = j, l1 = NA)
      }
      
      other_df <- other_df |> add_row(method = "lp_tilted", init = i, time = NA, mem = NA)
    }
  }
  
  save(l1_df, other_df, file = "results/lp_tilted.RData")
} else if (method == "pep") {
  # Power EP
  source("scripts/pep.R")
  
  for (i in 1:n_init) {
    print(paste0("Testing initialisation ", i))
    
    r_init <- init_list[[i]]$r_init
    Q_init <- init_list[[i]]$Q_init
    res <- tryCatch({set.seed(1); pep(X, y, mu_beta, Sigma_beta, 
                                      alpha, r_init, Q_init, 
                                      min_pass, max_pass, thresh, verbose, F)},
                    error = function(e) NA)
    if (is.list(res)) {
      for (j in 1:p) {
        values <- dnorm(grid_points[j, ], res$mu[j], sqrt(res$Sigma[j, j]))
        l1 <- 1 - 0.5*trapz(grid_points[j, ], abs(mcmc_values[j, ] - values))
        l1_df <- l1_df |> add_row(method = "pep", init = i, j = j, l1 = l1)
      }
      
      mark_res_1 <- mark({set.seed(1); pep(X, y, mu_beta, Sigma_beta, 
                                           alpha, r_init, Q_init, 
                                           min_pass, max_pass, thresh, verbose, F)})
      
      mark_res_2 <- mark({set.seed(1); pep(X, y, mu_beta, Sigma_beta, 
                                           alpha, r_init, Q_init, 
                                           min_pass, 1, thresh, verbose, T)})
      
      time <- as.numeric(mark_res_1$median)
      mem <- as.numeric(mark_res_2$mem_alloc)
      
      other_df <- other_df |> add_row(method = "pep", init = i, time = time, mem = mem)
    } else {
      for (j in 1:p) {
        l1_df <- l1_df |> add_row(method = "pep", init = i, j = j, l1 = NA)
      }
      
      other_df <- other_df |> add_row(method = "pep", init = i, time = NA, mem = NA)
    }
  }
  
  save(l1_df, other_df, file = "results/pep.RData")
} else if (method == "dep") {
  # Damped EP
  source("scripts/dep.R")
  
  for (i in 1:n_init) {
    print(paste0("Testing initialisation ", i))
    
    r_init <- init_list[[i]]$r_init
    Q_init <- init_list[[i]]$Q_init
    res <- tryCatch({set.seed(1); dep(X, y, mu_beta, Sigma_beta, 
                                      beta, r_init, Q_init, 
                                      min_pass, max_pass, thresh, verbose, F)},
                    error = function(e) NA)
    if (is.list(res)) {
      for (j in 1:p) {
        values <- dnorm(grid_points[j, ], res$mu[j], sqrt(res$Sigma[j, j]))
        l1 <- 1 - 0.5*trapz(grid_points[j, ], abs(mcmc_values[j, ] - values))
        l1_df <- l1_df |> add_row(method = "dep", init = i, j = j, l1 = l1)
      }
      
      mark_res_1 <- mark({set.seed(1); dep(X, y, mu_beta, Sigma_beta, 
                                           beta, r_init, Q_init, 
                                           min_pass, max_pass, thresh, verbose, F)})
      
      mark_res_2 <- mark({set.seed(1); dep(X, y, mu_beta, Sigma_beta, 
                                           beta, r_init, Q_init, 
                                           min_pass, 1, thresh, verbose, T)})
      
      time <- as.numeric(mark_res_1$median)
      mem <- as.numeric(mark_res_2$mem_alloc)
      
      other_df <- other_df |> add_row(method = "dep", init = i, time = time, mem = mem)
    } else {
      for (j in 1:p) {
        l1_df <- l1_df |> add_row(method = "dep", init = i, j = j, l1 = NA)
      }
      
      other_df <- other_df |> add_row(method = "dep", init = i, time = NA, mem = NA)
    }
  }
  
  save(l1_df, other_df, file = "results/dep.RData")
} else if (method == "dlep") {
  # Double-loop EP
  source("scripts/dlep.R")
  
  for (i in 1:n_init) {
    print(paste0("Testing initialisation ", i))
    
    r_init <- init_list[[i]]$r_init
    Q_init <- init_list[[i]]$Q_init
    res <- tryCatch({set.seed(1); dlep(X, y, mu_beta, Sigma_beta, 
                                       min_inner, max_inner, thresh_grad, outer_freq, r_init, Q_init, 
                                       min_pass, max_pass_dl, thresh, verbose, F)},
                    error = function(e) NA)
    if (is.list(res)) {
      for (j in 1:p) {
        values <- dnorm(grid_points[j, ], res$mu[j], sqrt(res$Sigma[j, j]))
        l1 <- 1 - 0.5*trapz(grid_points[j, ], abs(mcmc_values[j, ] - values))
        l1_df <- l1_df |> add_row(method = "dlep", init = i, j = j, l1 = l1)
      }
      
      mark_res_1 <- mark({set.seed(1); dlep(X, y, mu_beta, Sigma_beta, 
                                            min_inner, max_inner, thresh_grad, outer_freq, r_init, Q_init, 
                                            min_pass, max_pass_dl, thresh, verbose, F)})
      
      mark_res_2 <- mark({set.seed(1); dlep(X, y, mu_beta, Sigma_beta, 
                                            min_inner, max_inner, thresh_grad, outer_freq, r_init, Q_init, 
                                            min_pass, 1, thresh, verbose, T)})
      
      time <- as.numeric(mark_res_1$median)
      mem <- as.numeric(mark_res_2$mem_alloc)
      
      other_df <- other_df |> add_row(method = "dlep", init = i, time = time, mem = mem)
    } else {
      for (j in 1:p) {
        l1_df <- l1_df |> add_row(method = "dlep", init = i, j = j, l1 = NA)
      }
      
      other_df <- other_df |> add_row(method = "dlep", init = i, time = NA, mem = NA)
    }
  }
  
  save(l1_df, other_df, file = "results/dlep.RData")
} else if (method == "sep") {
  # Stochastic EP
  source("scripts/sep.R")
  
  for (i in 1:n_init) {
    print(paste0("Testing initialisation ", i))
    
    r_init <- init_list[[i]]$r_init
    Q_init <- init_list[[i]]$Q_init
    res <- tryCatch({set.seed(1); sep(X, y, mu_beta, Sigma_beta, 
                                      r_init, Q_init, 
                                      min_pass, max_pass, thresh, verbose, F)},
                    error = function(e) NA)
    if (is.list(res)) {
      for (j in 1:p) {
        values <- dnorm(grid_points[j, ], res$mu[j], sqrt(res$Sigma[j, j]))
        l1 <- 1 - 0.5*trapz(grid_points[j, ], abs(mcmc_values[j, ] - values))
        l1_df <- l1_df |> add_row(method = "sep", init = i, j = j, l1 = l1)
      }
      
      mark_res_1 <- mark({set.seed(1); sep(X, y, mu_beta, Sigma_beta, 
                                           r_init, Q_init, 
                                           min_pass, max_pass, thresh, verbose, F)})
      
      mark_res_2 <- mark({set.seed(1); sep(X, y, mu_beta, Sigma_beta, 
                                           r_init, Q_init, 
                                           min_pass, 1, thresh, verbose, T)})
      
      time <- as.numeric(mark_res_1$median)
      mem <- as.numeric(mark_res_2$mem_alloc)
      
      other_df <- other_df |> add_row(method = "sep", init = i, time = time, mem = mem)
    } else {
      for (j in 1:p) {
        l1_df <- l1_df |> add_row(method = "sep", init = i, j = j, l1 = NA)
      }
      
      other_df <- other_df |> add_row(method = "sep", init = i, time = NA, mem = NA)
    }
  }
  
  save(l1_df, other_df, file = "results/sep.RData")
} else if (method == "bbep") {
  # Black-box EP
  source("scripts/bbep.R")
  
  for (i in 1:n_init) {
    print(paste0("Testing initialisation ", i))
    
    r_init <- init_list[[i]]$r_init
    Q_init <- init_list[[i]]$Q_init
    res <- tryCatch({set.seed(1); bbep(X, y, mu_beta, Sigma_beta, 
                                       M_bb, epsilon, tau, r_init, Q_init, 
                                       min_iter, max_iter, thresh_grad, verbose, F)},
                    error = function(e) NA)
    if (is.list(res)) {
      for (j in 1:p) {
        values <- dnorm(grid_points[j, ], res$mu[j], sqrt(res$Sigma[j, j]))
        l1 <- 1 - 0.5*trapz(grid_points[j, ], abs(mcmc_values[j, ] - values))
        l1_df <- l1_df |> add_row(method = "bbep", init = i, j = j, l1 = l1)
      }
      
      mark_res_1 <- mark({set.seed(1); bbep(X, y, mu_beta, Sigma_beta, 
                                            M_bb, epsilon, tau, r_init, Q_init, 
                                            min_iter, max_iter, thresh_grad, verbose, F)})
      
      mark_res_2 <- mark({set.seed(1); bbep(X, y, mu_beta, Sigma_beta, 
                                            M_bb, epsilon, tau, r_init, Q_init, 
                                            min_iter, 1, thresh_grad, verbose, T)})
      
      time <- as.numeric(mark_res_1$median)
      mem <- as.numeric(mark_res_2$mem_alloc)
      
      other_df <- other_df |> add_row(method = "bbep", init = i, time = time, mem = mem)
    } else {
      for (j in 1:p) {
        l1_df <- l1_df |> add_row(method = "bbep", init = i, j = j, l1 = NA)
      }
      
      other_df <- other_df |> add_row(method = "bbep", init = i, time = NA, mem = NA)
    }
  }
  
  save(l1_df, other_df, file = "results/bbep.RData")
} else if (method == "auto_ep") {
  # Automatic tilted inference
  source("scripts/auto_ep.R")
  
  set.seed(1)
  train_inds <- sample(1:nrow(X), 0.8*nrow(X))
  train_df <- gen_train(X[train_inds, ], y[train_inds], mu_beta, Sigma_beta, 
                        M_auto, r_init, Q_init, 
                        min_pass, max_pass, thresh, Inf, verbose)
  
  for (i in 1:n_init) {
    print(paste0("Testing initialisation ", i))
    
    r_init <- init_list[[i]]$r_init
    Q_init <- init_list[[i]]$Q_init
    
    u_max <- knn_u_max(train_df[, c("mu", "sigma_2")], train_df$mu_h, train_df$sigma_2_h, D_max, D_max, k, folds)
    
    res <- tryCatch({set.seed(1); auto_ep(X, y, mu_beta, Sigma_beta, 
                                          k, train_df, u_max, M_auto, r_init, Q_init, 
                                          min_pass, max_pass, thresh, Inf, verbose, F)},
                    error = function(e) NA)
    if (is.list(res)) {
      for (j in 1:p) {
        values <- dnorm(grid_points[j, ], res$mu[j], sqrt(res$Sigma[j, j]))
        l1 <- 1 - 0.5*trapz(grid_points[j, ], abs(mcmc_values[j, ] - values))
        l1_df <- l1_df |> add_row(method = "auto_ep", init = i, j = j, l1 = l1)
      }
      
      mark_res_1 <- mark({set.seed(1); auto_ep(X, y, mu_beta, Sigma_beta, 
                                               k, train_df, u_max, M_auto, r_init, Q_init, 
                                               min_pass, max_pass, thresh, Inf, verbose, F)})
      
      mark_res_2 <- mark({set.seed(1); auto_ep(X, y, mu_beta, Sigma_beta, 
                                               k, train_df, u_max, M_auto, r_init, Q_init, 
                                               min_pass, 1, thresh, Inf, verbose, T)})
      
      time <- as.numeric(mark_res_1$median)
      mem <- as.numeric(mark_res_2$mem_alloc)
      
      other_df <- other_df |> add_row(method = "auto_ep", init = i, time = time, mem = mem)
    } else {
      for (j in 1:p) {
        l1_df <- l1_df |> add_row(method = "auto_ep", init = i, j = j, l1 = NA)
      }
      
      other_df <- other_df |> add_row(method = "auto_ep", init = i, time = NA, mem = NA)
    }
  }
  
  save(l1_df, other_df, file = "results/auto_ep.RData")
} else {
  stop("Invalid method")
}
