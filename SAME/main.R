set.seed(2020)
set.seed(2021)
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
same_input <- generate_same_input(T,D,K,pi_ber,N,Iteration,corrupt_pi=0.2)
########

Y0 = as.matrix(same_input$Y0)
X = same_input$X
T = same_input$T
K = same_input$K
D = same_input$D
c_k = same_input$c_k
C0 = sum(c_k)
W_tilde = same_input$W_tilde # observed w_tilde
YSG =  same_input$YSG
true_z = same_input$true_Z
true_w =  same_input$true_w$w
raw_X = same_input$raw_X # observation may after corruption
true_v = same_input$true_w$v
true_gamma = same_input$true_w$gamma
original_X = same_input$X_original # single cell before corruption

## check out assumption
#true_w = 0.5*W_tilde + 0.5*true_w[[1]]
#true_e  = Y0-true_w%*%true_z
#str(true_e)
#hist(true_e,100)
#hist(Y0,100)

cbind(raw_X[[1]]$w_tilde[,2],W_tilde[,2])
cbind(original_X[[1]]$w_tilde[,2],W_tilde[,2],raw_X[[1]]$w_tilde[,2])

# Starting values
mcmc_samples_theta1 = 50
Lambda = c(0:mcmc_samples_theta1) # Lambda = c(0,1,2,3,...,100)
#Lambda = c(0,rep(1,mcmc_samples_theta1))


rst <- SAME(Y0, X, W_tilde,
                 mcmc_samples_theta1, Lambda, c_k, YSG, alpha = 1)


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

i=1
cbind(gamma_est[[20]][,i]*v_est[[20]][,i],true_gamma[,i]*true_v[,i])
tbl <- table(gamma_est[[20]][,i],true_gamma[,i])
print(tbl)
chisq.test(tbl)
plot(v_est[[20]][which(true_gamma[,i]==1),i],true_v[which(true_gamma[,i]==1),i], xlab = 'estmation',ylab='true')
plot(v_est[[20]][which(gamma_est[[20]][,i]==1),i],true_v[which(gamma_est[[20]][,i]==1),i], xlab = 'estmation',ylab='true', main='v')
cor(v_est[[20]][which(gamma_est[[20]][,i]==1),i],true_v[which(gamma_est[[20]][,i]==1),i])

cbind(v_est[[20]][,i],true_v[,i])
i=1
plot(w_est[[i]],true_w[[i]], xlab = 'estmation',ylab='true', main="W t=1")
i=2
plot(w_est[[i]],true_w[[i]], xlab = 'estmation',ylab='true', main="W t=2")

plot(z_est,true_z, xlab = 'estmation',ylab='true', main = 'Z')
boxplot(as.vector(z_est-true_z))

ratio <- true_w[[i]] / w_est[[i]]
hist(ratio,100)$mids[c(61,67,86)]
## 1.21 1.33 1.71
## 1.2 - 1.22
## 1.32 - 1.34
## 1.70 - 1.72
idx = which(ratio >1.2, arr.ind = T)
new_idx = matrix(0,nrow = 1, ncol=2)
for(s in 1:nrow(idx)){
    if(ratio[idx[s,1],idx[s,2]] < 1.22) new_idx <- rbind(new_idx, idx[s,])
}
new_idx <- new_idx[-1,]


## T = 50
saveRDS(same_input, "write/simulation/same_input_T.8.v2020.rds")
saveRDS(rst, "write/simulation/rst_T.8.v2020.rds")
#rst.alpha05 <- readRDS("write/simulation/rst_T.15.v2020.rds" )
z_est_alpha05 <- rst$theta1$z
z_est_alpha1 <- rst$theta1$z

compare = data.frame(diff = c(as.vector(z_est_alpha1-true_z), as.vector(z_est_alpha05-true_z)), 
                method = c(rep("alpha 1",K*N),rep("alpha 0.5",K*N)))
boxplot(diff~method, data=compare)