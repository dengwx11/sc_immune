options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

T <- as.numeric(args[1])
D <- as.numeric(args[2])
K <- as.numeric(args[3])
alpha <- as.numeric(args[4])
corrupt_pi <- as.numeric(args[5])
tau_w_para <- as.numeric(args[6])
tau_xd_beta_para <- as.numeric(args[7])


set.seed(2021)
library(nnls)
library(MuSiC)
library(xbioc)
library(ggplot2)
library(Biobase)
library(plotROC)
#set.seed(2021)
source('SAME/run_same.R')
source('SAME/run_Simulation.R')

######### semi simulation
# same_input <- readRDS('data/Pseudo_Bulk/SAME_Input.rds')
################### simulation
######## tune parameter
# T=5
# D=500
# K=8
pi_ber = 0.3
N = 200 # bulk Y sample size
Iteration = 500 ## iteration number to get the largest angle between the vectors

# tau_w_para = 1  ## var(w)= 1/tau_w_para, smaller value --> bigger tissue-specific differentiation
# tau_xd_beta_para = 1  ## var(x)=tau_xd_beta_para approximately, 
#                        ## tau_xd = rgamma(D, 1, tau_xd_beta_para)
#                        ## smaller value --> smaller single cell expression variation
# corrupt_pi =  0.8 ## corruption rate due to shallow sequencing depth
########
set.seed(2021)
same_input <- generate_same_input(T,D,K,pi_ber,N,Iteration,corrupt_pi=corrupt_pi)
str_para = paste0("T=",T,".D=",D,".K=",K,".corrupt=",corrupt_pi,".tauW=",tau_w_para,".tauXdBeta=",tau_xd_beta_para)


Y0 = as.matrix(same_input$Y0)
X = same_input$X  ## seurat type for observed expression profile
T = same_input$T
K = same_input$K
D = same_input$D
c_k = same_input$c_k
C0 = sum(c_k)
W_tilde = same_input$W_tilde ## observed w_tilde
YSG =  same_input$YSG
true_z = same_input$true_Z
true_w =  same_input$true_w$w ## true w for T tissues
raw_X = same_input$raw_X ## list type for observation after corruption 
true_v = same_input$true_w$v 
true_gamma = same_input$true_w$gamma
# original_X = same_input$X_original # single cell before corruption

cbind(true_v*true_gamma[,1],true_w[[1]][,1],raw_X[[1]]$w_tilde[,1],true_w[[2]][,1],raw_X[[2]]$w_tilde[,1],W_tilde[,1])

# Starting values
mcmc_samples_theta1 = 100
Lambda = c(0:mcmc_samples_theta1) # Lambda = c(0,1,2,3,...,100)

alpha=1
rst1 <- SAME(Y0, X, W_tilde,
            mcmc_samples_theta1, Lambda, c_k, YSG, alpha =1)
rst$same_rst <- rst1
rst$same_input <- same_input
mydir <- "/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/simulation_rst/"            
saveRDS(rst,paste0(mydir,"rst.alpha=",alpha,'.',str_para,'.rds'))
