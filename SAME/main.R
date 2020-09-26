set.seed(2020)
source('run_same.R')


## Input (upper case)
## Y0: bulk for tissue t0, D*N
## X = list() : single cell for all tissues, tissue t0 is the first item
## W_tilde : empirical signature matrix for tissue t0, D*K
same_input <- readRDS('data/Pseudo_Bulk/SAME_Input.rds')
Y0 = as.matrix(same_input$Y0)
X = same_input$X
T = same_input$T
K = same_input$K
D = same_input$D
c_k = same_input$c_k
C0 = sum(c_k)
W_tilde = same_input$W_tilde
#YSG =  same_input$YSG


# Starting values
mcmc_samples_theta1 = 100
Lambda = c(0:mcmc_samples_theta1)

rst <- SAME(Y0, X, W_tilde,
                 mcmc_samples_theta1, Lambda)
