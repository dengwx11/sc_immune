set.seed(2020)
source("Update_theta1.R")
source("Update_theta2.R")

## Input (upper case)
## Y0: bulk for tissue t0
## X = list() : single cell for all tissues, tissue t0 is the first item
## W_tilde : empirical signature matrix for tissue t0

## Initialization (lower case)
# Noninformative Prior Parameters
tau_v = 0.01
alpha_prior_e = 10^(-6)
beta_prior_e = 10^(-6)
alpha_prior_x = 10^(-6)
beta_prior_x = 10^(-6)
alpha_prior_w = 10^(-6)
beta_prior_w = 10^(-6)
alpha_prior_pi = 10^(-6)
beta_prior_pi = 10^(-6)

# Starting values
mcmc_samples_theta2 = 100
Lambda = c(0:mcmc_samples_theta2)
D = 100 # number of genes
K = 20 # number of cell types
T = 30  # number of tissues
N = 100 # number of bulk samples
c_k = matrix(200,nrow = K, ncol = T) # sequence of cell counts for each cell type and each tissue
C0 = sum(c_k) # total number of single cell counts

# Estimations
# numerber or vector is stored as mcmc sequence
# matrix is stored as the newest one 
# theta2
mcmc_samples_theta1 = sum(Lambda)
tau_e_est <- matrix(0, nrow =  mcmc_samples_theta1, ncol = D)
alpha_unif_est <- matrix(0, nrow =  mcmc_samples_theta1, ncol = 1)
gamma_est <- matrix(rbinom(D*K,1,0.5),nrow = D, ncol = K)
v_est <- matrix(rnorm(D*K,mean = 2,sd = 1),nrow = D, ncol = K)
pi_ber_est <- matrix(0, nrow =  mcmc_samples_theta1, ncol = 1)
tau_x_est <- matrix(0, nrow =  mcmc_samples_theta1, ncol = D)
tau_w_est <- matrix(0, nrow =  mcmc_samples_theta1, ncol = 1)

# theta1
z_est <- matrix(rnorm(K*N,mean = 2,sd = 1),nrow = K, ncol = N)
w_est <- list()
for(i in 1:T) w_est[[i]] <- matrix(rbinom(D*K,1,0.5),nrow = D, ncol = K) # first item = t0 tissue 


## SAME

for(i in 1:mcmc_samples){
    ## Step 1: Update theta2 for Lambda[i] times
    start_idx = sum(Lambda[1:i])
    end_idx = sum(Lambda[1:(i+1)])
    for(j in start_idx:end_idx){

    }

    ## Step 2: Update theta1 for one time

}

