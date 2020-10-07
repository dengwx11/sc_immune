set.seed(2020)
source('SAME/run_same.R')
source('SAME/run_Simulation.R')


## Input (upper case)
## Y0: bulk for tissue t0, D*N
## X = list() : single cell for all tissues, tissue t0 is the first item
## W_tilde : empirical signature matrix for tissue t0, D*K

######### semi simulation
same_input <- readRDS('data/Pseudo_Bulk/SAME_Input.rds')
######## simulation
T=4
D=500
K=3
pi_ber = 0.3
N = 200 # bulk Y sample size
Iteration = 500 ## iteration number to get the largest angle between the vectors
same_input <- generate_same_input(T,D,K,pi_ber,N,Iteration)
########

Y0 = as.matrix(same_input$Y0)
X = same_input$X
T = same_input$T
K = same_input$K
D = same_input$D
c_k = same_input$c_k
C0 = sum(c_k)
W_tilde = same_input$W_tilde
YSG =  same_input$YSG
true_z = same_input$true_Z
true_w_sc =  same_input$true_w$w
raw_X = same_input$raw_X


## check out assumption
true_w = 0.5*W_tilde + 0.5*true_w_sc[[1]]
true_e  = Y0-true_w%*%true_z
str(true_e)
hist(true_e,100)
hist(Y0,100)

cbind(raw_X[[1]]$w_tilde[,2],W_tilde[,2])

# Starting values
mcmc_samples_theta1 = 100
Lambda = c(0:mcmc_samples_theta1) # Lambda = c(0,1,2,3,...,100)
Lambda = c(0,rep(1,mcmc_samples_theta1))


rst <- SAME(Y0, X, W_tilde,
                 mcmc_samples_theta1, Lambda, c_k, YSG)
