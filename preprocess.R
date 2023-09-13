
# Data preprocessing

source("setup.R")
set.seed(1)

## Generating the data

o_rings <- read.csv("data/challenger.csv")
X <- unname(cbind(1, scale(o_rings[, 2:3])))
y <- as.numeric(o_rings[, 1] != 0)

N <- nrow(X)
p <- ncol(X)

mu_beta <- rep(0, p)
Sigma_beta <- 10000*diag(p)

save(N, p, X, y, mu_beta, Sigma_beta, file = "data/o_rings.RData")

## MCMC gold standard

stan_res <- stan(file = "scripts/mcmc.stan", 
                 data = list(N = N, p = p, X = X, y = y,
                             mu_beta = mu_beta,
                             Sigma_beta = Sigma_beta),
                 chains = mcmc_chains,
                 cores = mcmc_cores,
                 warmup = mcmc_warmup,
                 iter = mcmc_iter, 
                 refresh = 100,
                 init = 0)

mcmc_samples <- extract(stan_res)$beta
mcmc_mu <- colMeans(mcmc_samples)
mcmc_Sigma <- var(mcmc_samples)
mcmc_rhat_max <- max(summary(stan_res)$summary[paste0("beta[", 1:p, "]"), "Rhat"])

if (mcmc_rhat_max > rhat_thresh) {
  stop("R-hat is too high")
}

grid_points <- matrix(nrow = p, ncol = total_grid_points)
mcmc_values <- matrix(nrow = p, ncol = total_grid_points)

for (j in 1:p) {
  density_res <- density(mcmc_samples[, j], bw = "SJ-ste",
                         from = mcmc_mu[j] - sd_multiple*sqrt(mcmc_Sigma[j, j]),
                         to = mcmc_mu[j] + sd_multiple*sqrt(mcmc_Sigma[j, j]),
                         n = total_grid_points)
  
  grid_points[j, ] <- density_res$x
  mcmc_values[j, ] <- density_res$y
}

save(grid_points, mcmc_values, file = "results/mcmc.RData")
