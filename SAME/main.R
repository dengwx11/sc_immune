set.seed(2020)
source('run_same.R')


## Input (upper case)
## Y0: bulk for tissue t0, D*N
## X = list() : single cell for all tissues, tissue t0 is the first item
## W_tilde : empirical signature matrix for tissue t0, D*K

# Starting values
mcmc_samples_theta1 = 100
Lambda = c(0:mcmc_samples_theta1)

rst <- SAME(Y0, X, W_tilde,
                 mcmc_samples_theta1, Lambda)
