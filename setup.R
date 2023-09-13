
# Set-up

library(tidyverse)
library(rstan)
library(kernlab)
library(bench)
library(doParallel)
library(abind)
library(sn)
library(plyr)
library(condmixt)
library(Matrix)
library(calculus)
library(pracma)
library(matrixsampling)

mcmc_chains <- 4
mcmc_cores <- 4
mcmc_warmup <- 1000
mcmc_iter <- 10000
rhat_thresh <- 1.1

total_grid_points <- 1000
sd_multiple <- 5
n_init <- 10

M_bb <- 50
M_auto <- 1000
epsilon <- 0.01
tau <- 50
alpha <- 0.5
beta <- 0.5
min_pass <- 1
min_inner <- 10
min_iter <- 1
max_pass <- 50
max_pass_dl <- 500
max_inner <- 200
max_iter <- 50
outer_freq <- 1
k <- 5
folds <- 10
D_max <- 0.01
thresh <- 0.01
thresh_grad <- 0.001
verbose <- F
