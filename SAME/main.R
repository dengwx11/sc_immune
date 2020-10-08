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
T=100
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
true_w =  same_input$true_w$w
raw_X = same_input$raw_X
true_v = same_input$true_w$v
true_gamma = same_input$true_w$gamma

## check out assumption
true_w = 0.5*W_tilde + 0.5*true_w[[1]]
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


z_est <-   rst$theta1$z 
w_est <-  rst$theta1$w  
tau_e_est <-  rst$theta2$tau_e  
alpha_unif_est  <- rst$theta2$alpha_unif  
gamma_est  <-  rst$theta2$gamma  
v_est  <-   rst$theta2$v 
pi_ber_est  <-  rst$theta2$pi_ber  
tau_x_est  <-  rst$theta2$tau_x  
tau_w_est  <-   rst$theta2$tau_w 


## for debug
#hist((v_est[[k]]*gamma_est[[k]])[which(gamma_est[[k]]==1)],100)

# i=2
# cbind(gamma_est[[20]][,i]*v_est[[20]][,i],true_gamma[,i]*true_v[,i])
# table(gamma_est[[20]][,i],true_gamma[,i])
# plot(v_est[[20]][which(true_gamma[,i]==1),i],true_v[which(true_gamma[,i]==1),i], xlab = 'estmation',ylab='true')
# cbind(v_est[[20]][,i],true_v[,i])
# i=3
# plot(w_est[[i]],true_w[[i]], xlab = 'estmation',ylab='true')

# plot(z_est,true_z, xlab = 'estmation',ylab='true')