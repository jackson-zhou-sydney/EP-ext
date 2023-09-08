
# Set-up

library(tidyverse)
library(rstan)
library(kernlab)
library(doParallel)
library(abind)
library(sn)

registerDoParallel(cores = 8)

mcmc_chains <- 4
mcmc_cores <- 4
mcmc_warmup <- 1000
mcmc_iter <- 10000
rhat_thresh <- 1.1
total_grid_points <- 1000
sd_multiple <- 5
