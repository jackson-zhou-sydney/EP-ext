
# Numerical experiments

source("setup.R")
load("data/musk.RData")
load("results/mcmc.RData")
method <- "ep"

if (method == "ep") {
  # Standard EP
} else if (method == "parallel_ep") {
  # Distributed computation
} else if (method == "lp") {
  # Laplace propagation
} else if (method == "power_ep") {
  # Power EP
} else if (method == "damped_ep") {
  # Damped EP
} else if (method == "dl_ep") {
  # Double-loop EP
} else if (method == "sep") {
  # Stochastic EP
} else if (method == "bb_ep") {
  # Black-box EP
} else if (method == "auto_ep") {
  # Automatic tilted inference
} else {
  stop("Invalid method")
}
