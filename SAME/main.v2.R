set.seed(2021)
library(nnls)
library(MuSiC)
library(xbioc)
library(ggplot2)
library(Biobase)
library(plotROC)
#set.seed(2021)
source('SAME/run_same.v2.R')
source('SAME/run_Simulation.R')

######### semi simulation
# same_input <- readRDS('data/Pseudo_Bulk/SAME_Input.rds')
################### simulation
######## tune parameter
T=5
D=500
K=8
pi_ber = 0.3
N = 200 # bulk Y sample size
Iteration = 500 ## iteration number to get the largest angle between the vectors

tau_w_para = 1  ## var(w)= 1/tau_w_para, smaller value --> bigger tissue-specific differentiation
tau_xd_beta_para = 1  ## var(x)=tau_xd_beta_para approximately, 
                       ## tau_xd = rgamma(D, 1, tau_xd_beta_para)
                       ## smaller value --> smaller single cell expression variation
corrupt_pi =  0.4 ## corruption rate due to shallow sequencing depth
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
mcmc_samples_theta1 = 30
Lambda = c(0:mcmc_samples_theta1) # Lambda = c(0,1,2,3,...,100)

alpha=1
empirical_pi=0.3
rst1 <- SAME(X, empirical_pi,
            mcmc_samples_theta1, Lambda, c_k, YSG, alpha =1)


w_est_alpha1 <-  rst1$theta1$w
est_vg <- lapply(c(1:mcmc_samples_theta1), function(i) rst1$theta2$gamma[[i]] * rst1$theta2$v[[i]])
average_est_vg_1 <- Reduce('+',est_vg)/mcmc_samples_theta1            

true_vg <- true_v * true_gamma

i=1
plot(w_est_alpha1[[i]],true_vg, xlab = 'estmation',ylab='true', main="w_est_alpha1 t=1")
plot(average_est_vg_1,true_vg, xlab = 'estmation',ylab='true', main="v*g t=1")
plot(W_tilde,true_vg, xlab = 'W tilde',ylab='true', main="W_tilde t=1")
plot(W_tilde,w_est_alpha1[[1]], xlab = 'W tilde',ylab='est', main="W_tilde t=1")

compare.w = data.frame(diff = c(as.vector(w_est_alpha1[[i]]-true_vg),
                                as.vector(average_est_vg_1-true_vg),
                                as.vector(W_tilde-true_vg)), 
                     method = c(rep("w_est",K*D),rep("est_vg",K*D),rep("W_tilde",K*D)))
#boxplot(diff~method, data=compare.w)
ggplot(compare.w, aes(x=method,y=diff)) +  geom_boxplot()+ ylim(-1,2) 

## evaluate z
z_w_tilde_nnls  <- sapply(c(1:N), function(j) nnls(W_tilde, Y0[,j])$x)
z_est_vg_nnls  <- sapply(c(1:N), function(j) nnls(average_est_vg_1, Y0[,j])$x)

plot(z_w_tilde_nnls,true_z, xlab = 'estmation',ylab='true', main = 'Z')
plot(z_est_vg_nnls,true_z, xlab = 'estmation',ylab='true', main = 'Z')


##### music ####
assaycounts <- raw_X[[1]]$counts
rownames(assaycounts) <- paste0('Gene',c(1:D))
colnames(assaycounts) <- paste0('Cell',c(1:ncol(assaycounts)))
pheno <- X[[1]]@meta.data
pheno$SampleID <- paste0('Cell',c(1:ncol(assaycounts)))
rownames(pheno) <- paste0('Cell',c(1:ncol(assaycounts)))
sc.eset <- ExpressionSet(assayData = assaycounts, phenoData = as(data.frame(pheno),"AnnotatedDataFrame"))
assaycounts <- Y0
rownames(assaycounts) <- paste0('Gene',c(1:D))
colnames(assaycounts) <- paste0('ID',c(1:N))
pheno <- data.frame('subjectID' = paste0('ID',c(1:N)))
rownames(pheno) <- paste0('ID',c(1:N))
bulk.eset <- ExpressionSet(assayData = assaycounts, phenoData = as(data.frame(pheno),"AnnotatedDataFrame"))
# Estimate cell type proportions
Est.prop = music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'Celltype_used',
                               samples = 'SampleID',  verbose = F)
names(Est.prop)
z_est_music <- t(Est.prop$Est.prop.weighted)
z_est_nnls <- t(Est.prop$Est.prop.allgene)

cor(as.vector(z_est_vg_nnls),as.vector(true_z))
cor(as.vector(z_w_tilde_nnls),as.vector(true_z))
cor(as.vector(z_est_music),as.vector(true_z))
cor(as.vector(z_est_nnls),as.vector(true_z))

compare = data.frame(diff = c(as.vector(z_est_vg_nnls-true_z),as.vector(z_w_tilde_nnls-true_z), as.vector(z_est_nnls-true_z),
                              as.vector(z_est_music-true_z)), 
                     method = c(rep("est_vg",K*N),rep("w_tilde_nnls",K*N),rep("MuSiC_NNLS",K*N),rep("MuSiC",K*N)),
                     para = "simulation11")   
ggplot(compare, aes(x=method,y=diff)) +  geom_boxplot()+ ylim(-.3,.4) 

## evaluate z
plot(z_est_vg_nnls,true_z, xlab = 'estmation',ylab='true', main = 'Z')
plot(z_est_music,true_z, xlab = 'estmation',ylab='true', main = 'Z')
plot(z_est_nnls,true_z, xlab = 'estmation',ylab='true', main = 'Z')


average_gamma_est_1 <- Reduce('+',rst1$theta2$gamma)/mcmc_samples_theta1
tbl_1 <- data.frame('Estmated_Gamma' = average_gamma_est_1[,i],"True_Gamma" = true_gamma[,i])
tbl_1$corruption_prob <- factor(corrupt_pi)
tbl_1$beta <- tau_xd_beta_para
tbl_1$Estimated_Gamma_bin <- 1*(tbl_1$Estmated_Gamma >= .5)
tbl_1$alpha <- 1
tbl_1$corrupt_pi <- corrupt_pi

tbl_1_02 <- tbl_1
tbl_1_04 <- tbl_1
tbl_1_06 <- tbl_1


ggplot(tbl_1, aes(d=True_Gamma,m=Estmated_Gamma,color=factor(corrupt_pi)))+geom_roc(n.cuts = 0)

## EB
gamma <- get_est_vg(rst1)$averg_gamma
est_vg <- get_est_vg(rst1)$vg
# v
idc_W_tilde_orig <-(1*(est_vg>0))
sum(idc_W_tilde_orig)/sum(est_vg^2*idc_W_tilde_orig)
#pi
mean(gamma)
