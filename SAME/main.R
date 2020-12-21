options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

T <- as.numeric(args[1])
D <- as.numeric(args[2])
K <- as.numeric(args[3])
alpha <- as.numeric(args[4])




set.seed(2020)
library(nnls)
library(MuSiC)
library(xbioc)
library(ggplot2)
library(Biobase)
library(plotROC)
#set.seed(2021)
source('SAME/run_same.R')
source('SAME/run_Simulation.R')


## Input (upper case)
## Y0: bulk for tissue t0, D*N
## X = list() : single cell for all tissues, tissue t0 is the first item
## W_tilde : empirical signature matrix for tissue t0, D*K

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
tau_xd_beta_para = 1  ## var(x)=1/tau_xd_beta_para approximately, 
                       ## tau_xd = rgamma(D, 1, tau_xd_beta_para)
                       ## smaller value --> smaller single cell expression variation
corrupt_pi =  0.2 ## corruption rate due to shallow sequencing depth
########
set.seed(2020)
same_input <- generate_same_input(T,D,K,pi_ber,N,Iteration,corrupt_pi=corrupt_pi)
str_para = paste0("T=",T,".D=",D,".K=",K,".corrupt=",corrupt_pi,".tauW=",tau_w_para,".tauXdBeta=",tau_xd_beta_para)
###################



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
raw_X = same_input$raw_X ## list type for observation may after corruption 
true_v = same_input$true_w$v 
true_gamma = same_input$true_w$gamma
original_X = same_input$X_original # single cell before corruption

## check out assumption
#true_w = 0.5*W_tilde + 0.5*true_w[[1]]
#true_e  = Y0-true_w%*%true_z
#str(true_e)
#hist(true_e,100)
#hist(Y0,100)

# cbind(raw_X[[1]]$w_tilde[,2],W_tilde[,2])
cbind(original_X[[1]]$w_tilde[,1],W_tilde[,1],raw_X[[1]]$w_tilde[,1])

# Starting values
mcmc_samples_theta1 = 50
Lambda = c(0:mcmc_samples_theta1) # Lambda = c(0,1,2,3,...,100)
#Lambda = c(0,rep(1,mcmc_samples_theta1))



alpha=0.5
rst <- SAME(Y0, X, W_tilde,
            mcmc_samples_theta1, Lambda, c_k, YSG, alpha =alpha)
# saveRDS(rst,paste0('/Users/wenxuandeng/GoogleDrive/sucksalt/SC/sc_immune/code/sc_immune_copy/write/rst.alpha=',alpha,'.',str_para,'.rds'))
#average_gamma_est <- Reduce('+',gamma_est)/mcmc_samples_theta1


if(alpha==1){
  z_est_alpha1 <- rst$theta1$z
  z_est = z_est_alpha1
}else{
  z_est_alpha05 <- rst$theta1$z
  z_est = z_est_alpha05
}


alpha=.5
rst <- SAME(Y0, X, W_tilde,
            mcmc_samples_theta1, Lambda, c_k, YSG, alpha =alpha)
saveRDS(rst,paste0('/Users/wenxuandeng/GoogleDrive/sucksalt/SC/sc_immune/code/sc_immune_copy/write/rst.alpha=',alpha,'.',str_para,'.rds'))

if(alpha==1){
  z_est_alpha1 <- rst$theta1$z
  z_est = z_est_alpha1
}else{
  z_est_alpha05 <- rst$theta1$z
  z_est = z_est_alpha05
}


z_est <-   rst$theta1$z 
w_est <-  rst$theta1$w  
# tau_e_est <-  rst$theta2$tau_e  
# alpha_unif_est  <- rst$theta2$alpha_unif  
# gamma_est  <-  rst$theta2$gamma  
# v_est  <-   rst$theta2$v 
# pi_ber_est  <-  rst$theta2$pi_ber  
# tau_x_est  <-  rst$theta2$tau_x  
# tau_w_est  <-   rst$theta2$tau_w 
i=1
average_gamma_est <- Reduce('+',rst$theta2$gamma)/mcmc_samples_theta1
tbl0 <- data.frame('Estmated_Gamma' = average_gamma_est[,i],"True_Gamma" = true_gamma[,i])
tbl0$corruption_prob <- factor(corrupt_pi)
#tbl=tbl0
#tbl = rbind(tbl,tbl0)

# average_gamma_est <- Reduce('+',rst$theta2$gamma)/mcmc_samples_theta1
# tbl0 <- data.frame('Estmated_Gamma' = average_gamma_est[,i],"True_Gamma" = true_gamma[,i])
tbl0$beta <- tau_xd_beta_para

tbl=tbl0

# tbl <- rbind(tbl, tbl0)
# write.table(tbl,'./write/tbl.tauW=1.corrupt=0.2.tauXdBeta.txt',quote = F, row.names = F)

tbl = tbl[tbl$beta==0.5,]
tbl$tauW = 1
tbl0$tauW = 0.5
tbl <- rbind(tbl,tbl0)
write.table(tbl,'./write/tbl.tauXdBeta=0.5.tauW.txt',quote = F, row.names = F)
## for debug
#hist((v_est[[k]]*gamma_est[[k]])[which(gamma_est[[k]]==1)],100)

# i=1
# k=30
# average_gamma_est <- Reduce('+',gamma_est)/mcmc_samples_theta1
# cbind(gamma_est[[k]][,i]*v_est[[k]][,i],true_gamma[,i]*true_v[,i])
# tbl <- data.frame('Estmated_Gamma' = average_gamma_est[,i],"True_Gamma" = true_gamma[,i])
# boxplot(Estmated_Gamma~True_Gamma,data=tbl)
tbl$Estimated_Gamma_bin <- 1*(tbl$Estmated_Gamma >= .5)
# tbl$Method = 'alpha_1'
# print(tbl)
print(table(tbl$Estimated_Gamma_bin,tbl$True_Gamma))
# chisq.test(tbl)
# plot(v_est[[k]][which(true_gamma[,i]==1),i],true_v[which(true_gamma[,i]==1),i], xlab = 'estmation',ylab='true')
# plot(v_est[[k]][which(gamma_est[[k]][,i]==1),i],true_v[which(gamma_est[[k]][,i]==1),i], xlab = 'estmation',ylab='true', main='v')
# cor(v_est[[k]][which(gamma_est[[k]][,i]==1),i],true_v[which(gamma_est[[k]][,i]==1),i])

## W estimation evaluation
i=1
plot(w_est[[i]],true_w[[i]], xlab = 'estmation',ylab='true', main="W t=1")
i=2
plot(w_est[[i]],true_w[[i]], xlab = 'estmation',ylab='true', main="W t=2")

## ROC curve
# ROC.df <- data.frame('fdr' = 0,'tpr' = 0,"cutoff"=0)
# for(cutoff in c(1:100)/100){
#   Estimated_Gamma_bin <- 1*(tbl$Estmated_Gamma >= cutoff)
#   tbl.new <- table(Estimated_Gamma_bin,tbl$True_Gamma)
#   fdr <- tbl.new[2,1]/(D-sum(tbl$True_Gamma))
#   tpr <- tbl.new[2,2]/sum(tbl$True_Gamma)
#   ROC.df <- rbind(ROC.df, c(fdr,tpr,cutoff))
# }
# ROC.df <- ROC.df[-1,]
# plot(ROC.df$fdr,ROC.df$tpr)
# ROC.df$Method = "alpha_0.5"
# 
# ROC.df.1 <- data.frame('fdr' = 0,'tpr' = 0,"cutoff"=0)
# for(cutoff in c(1:100)/100){
#   Estimated_Gamma_bin <- 1*(tbl$Estmated_Gamma >= cutoff)
#   tbl.new <- table(Estimated_Gamma_bin,tbl$True_Gamma)
#   fdr <- tbl.new[2,1]/(D-sum(tbl$True_Gamma))
#   tpr <- tbl.new[2,2]/sum(tbl$True_Gamma)
#   ROC.df.1 <- rbind(ROC.df.1, c(fdr,tpr,cutoff))
# }
# ROC.df.1 <- ROC.df.1[-1,]
# plot(ROC.df.1$fdr,ROC.df.1$tpr)
# ROC.df.1$Method = "alpha_1"
# 
# ROC.df <- rbind(ROC.df,ROC.df.1)
ggplot(tbl, aes(d=True_Gamma,m=Estmated_Gamma,color=factor(beta)))+geom_roc(n.cuts = 0)
ggplot(tbl, aes(d=True_Gamma,m=Estmated_Gamma,color=factor(tauW)))+geom_roc(n.cuts = 0)
ggplot(tbl, aes(d=True_Gamma,m=Estmated_Gamma,color=corruption_prob))+labs(fill="corruption_prob")+geom_roc(n.cuts = 0)

# tau_xd
# plot(tau_x_est[nrow(tau_x_est),],raw_X[[1]]$tau_xd, xlab = 'estmation',ylab='true')
# plot(c(1:nrow(tau_w_est)),tau_w_est[,1])
# 
# 
# 
# cbind(v_est[[k]][,i],true_v[,i])
# i=1
# plot(w_est[[i]],true_w[[i]], xlab = 'estmation',ylab='true', main="W t=1")
# i=2
# plot(w_est[[i]],true_w[[i]], xlab = 'estmation',ylab='true', main="W t=2")
# 
plot(z_est,true_z, xlab = 'estmation',ylab='true', main = 'Z')
# boxplot(as.vector(z_est-true_z))
# 
# ratio <- true_w[[i]] / w_est[[i]]
# hist(ratio,100)$mids[c(61,67,86)]
## 1.21 1.33 1.71
## 1.2 - 1.22
## 1.32 - 1.34
## 1.70 - 1.72
# idx = which(ratio >1.2, arr.ind = T)
# new_idx = matrix(0,nrow = 1, ncol=2)
# for(s in 1:nrow(idx)){
#     if(ratio[idx[s,1],idx[s,2]] < 1.22) new_idx <- rbind(new_idx, idx[s,])
# }
# new_idx <- new_idx[-1,]


## T = 50
# saveRDS(same_input, "write/simulation/same_input_T.8.v2020.rds")
# saveRDS(rst, "write/simulation/rst_T.8.v2020.rds")
#rst.alpha05 <- readRDS("write/simulation/rst_T.15.v2020.rds" )



##### NNLS ####
 z_est_nnls  <- sapply(c(1:N), function(j) nnls(W_tilde, Y0[,j])$x)

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




##### CIBERSORTx ####
sigMat = data.frame("Gene symbol" = paste0('Gene',c(1:D)), W_tilde)
colnames(sigMat) <- c("Gene symbol",paste0('type',c(1:K)))
write.table(sigMat,file=paste0("./write/WTilde.",str_para,'.txt'),
            sep='\t',quote=F, row.names = F)
bulk = data.frame("Gene" = paste0('Gene',c(1:D)), Y0*20)
colnames(bulk) <- c("Gene",paste0('ID',c(1:N)))
write.table(bulk,file=paste0("./write/bulk.",str_para,'.txt'),
            sep="\t",quote=F, row.names = F)
sc = data.frame('Gene'=paste0('Gene',c(1:D)), raw_X[[1]]$counts)
colnames(sc) = c('Gene', raw_X[[1]]$Celltype_used)
write.table(sc,file=paste0("./write/sc.",str_para,'.txt'),
            sep='\t',quote=F, row.names = F)

z_est_cibersortx <-read.table('./write/CIBERSORTx_Job16_Results.txt', header=T,sep='\t')
z_est_cibersortx <- t(z_est_cibersortx[,c(2:(K+1))])/20

####### Final Comparison

compare = data.frame(diff = c(as.vector(z_est_alpha1-true_z), as.vector(z_est_alpha05-true_z), as.vector(z_est_nnls-true_z),
                              as.vector(z_est_music-true_z), as.vector(z_est_cibersortx-true_z)), 
                     method = c(rep("alpha_1",K*N),rep("alpha_0.5",K*N), rep("NNLS",K*N),rep("MuSiC",K*N),rep("CIBERSORTx",K*N)),
                     para = "simulation5")
compare = data.frame(diff = c(as.vector(z_est_alpha1-true_z), as.vector(z_est_alpha05-true_z), as.vector(z_est_nnls-true_z),
                              as.vector(z_est_music-true_z)), 
                     method = c(rep("alpha_1",K*N),rep("alpha_0.5",K*N), rep("NNLS",K*N),rep("MuSiC",K*N)),
                     para = "simulation11")
compare = data.frame(diff = c(as.vector(z_est_alpha05-true_z), as.vector(z_est_nnls-true_z),
                              as.vector(z_est_music-true_z)), 
                     method = c(rep("alpha_0.5",K*N), rep("NNLS",K*N),rep("MuSiC",K*N)),
                     para = "simulation11")                     
write.table(compare,'./write/simulation11.txt',quote = F,sep='\t',row.names = F)
## simulation1: T=5, D=500, K=8,mu=1,tau_v=4,tau_w=1,tau_xd =0.1
## simulation2: T=5, D=500, K=8,mu=1,tau_v=4,tau_w=1,tau_xd =0.5
## simulation3: T=5, D=500, K=8,mu=1,tau_v=4,tau_w=1,tau_xd =1
## simulation4: T=5, D=500, K=8,mu=1,tau_v=4,tau_w=.5,tau_xd =.5
## simulation5: T=5, D=500, K=8,mu=1,tau_v=4,tau_w=.5,tau_xd =.5, corruption = .2
## simulation6: T=5, D=500, K=8,mu=1,tau_v=4,tau_w=.5,tau_xd =1, corruption = .2
## simulation7: T=5, D=500, K=8,mu=1,tau_v=4,tau_w=.5,tau_xd =1, corruption = .4
## simulation8: T=5, D=500, K=8,mu=1,tau_v=4,tau_w=.5,tau_xd =.1, corruption = .4
## simulation9: T=5, D=500, K=8,mu=5,tau_v=1,tau_w=.5,tau_xd =.1, corruption = .4
## simulation10: T=5, D=500, K=8,mu=5,tau_v=1,tau_w=.5,tau_xd =.1, corruption = .6
## simulation11: T=5, D=500, K=8,mu=5,tau_v=1,tau_w=.5,tau_xd =.1, corruption = .8
boxplot(diff~method, data=compare)
median(z_est_alpha1-true_z)
median(z_est_alpha05-true_z)
median(z_est_nnls-true_z)
median(z_est_music-true_z)
median(z_est_cibersortx-true_z)

compare1 <- read.table('./write/simulation1.txt',sep='\t',header =T)
compare1$para <- "tau_w=1_tau_xd=0.1"
compare2 <- read.table('./write/simulation2.txt',sep='\t',header =T)
compare2$para <- "tau_w=1_tau_xd=0.5"
compare3 <- read.table('./write/simulation3.txt',sep='\t',header =T)
compare3$para <- "tau_w=1_tau_xd=1"
compare4 <- read.table('./write/simulation4.txt',sep='\t',header =T)
compare4$para <- "tau_w=.5_tau_xd=.5"
compare <- rbind(compare1,compare2,compare3,compare4)

ggplot(compare, aes(x=para,y=diff,fill=method))+geom_boxplot()
