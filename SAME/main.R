options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

T <- as.numeric(args[1])
D <- as.numeric(args[2])
K <- as.numeric(args[3])
alpha <- as.numeric(args[4])




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
tau_xd_beta_para = .1  ## var(x)=1/tau_xd_beta_para approximately, 
                       ## tau_xd = rgamma(D, 1, tau_xd_beta_para)
                       ## smaller value --> smaller single cell expression variation
corrupt_pi =  0.6 ## corruption rate due to shallow sequencing depth
########
set.seed(2021)
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
raw_X = same_input$raw_X ## list type for observation after corruption 
true_v = same_input$true_w$v 
true_gamma = same_input$true_w$gamma
original_X = same_input$X_original # single cell before corruption

cbind(true_w[[2]][,1],true_w[[1]][,1],original_X[[1]]$w_tilde[,1],W_tilde[,1])


## check out assumption
#true_w = 0.5*W_tilde + 0.5*true_w[[1]]
#true_e  = Y0-true_w%*%true_z
#str(true_e)
#hist(true_e,100)
#hist(Y0,100)

# cbind(raw_X[[1]]$w_tilde[,2],W_tilde[,2])


# Starting values
mcmc_samples_theta1 = 30
Lambda = c(0:mcmc_samples_theta1) # Lambda = c(0,1,2,3,...,100)
#Lambda = c(0,rep(1,mcmc_samples_theta1))



alpha=0.5
rst05 <- SAME(Y0, X, W_tilde,
            mcmc_samples_theta1, Lambda, c_k, YSG, alpha =.5)
# saveRDS(rst,paste0('/Users/wenxuandeng/GoogleDrive/sucksalt/SC/sc_immune/code/sc_immune_copy/write/rst.alpha=',alpha,'.',str_para,'.rds'))
#average_gamma_est <- Reduce('+',gamma_est)/mcmc_samples_theta1



z_est_alpha05 <- rst05$theta1$z
w_est_alpha05 <-  rst05$theta1$w



alpha=1
rst1 <- SAME(Y0, X, W_tilde,
            mcmc_samples_theta1, Lambda, c_k, YSG, alpha =1)
#saveRDS(rst,paste0('/Users/wenxuandeng/GoogleDrive/sucksalt/SC/sc_immune/code/sc_immune_copy/write/rst.alpha=',alpha,'.',str_para,'.rds'))


z_est_alpha1 <- rst1$theta1$z
w_est_alpha1 <-  rst1$theta1$w

alpha=0
rst0 <- SAME(Y0, X, W_tilde,
            mcmc_samples_theta1, Lambda, c_k, YSG, alpha =0)
#saveRDS(rst,paste0('/Users/wenxuandeng/GoogleDrive/sucksalt/SC/sc_immune/code/sc_immune_copy/write/rst.alpha=',alpha,'.',str_para,'.rds'))


z_est_alpha0 <- rst0$theta1$z
w_est_alpha0 <-  rst0$theta1$w

## evaluate z
plot(z_est_alpha05,true_z, xlab = 'estmation',ylab='true', main = 'Z')
plot(z_est_alpha1,true_z, xlab = 'estmation',ylab='true', main = 'Z')
plot(z_est_alpha0,true_z, xlab = 'estmation',ylab='true', main = 'Z')

compare.z = data.frame(diff = c(as.vector(z_est_alpha05-true_z), as.vector(z_est_alpha1-true_z),as.vector(z_est_alpha0-true_z)), 
                     method = c(rep("alpha_0.5",K*N),rep("alpha_1",K*N),rep("alpha_0",K*N)))
#boxplot(diff~method, data=compare.z) 
ggplot(compare.z, aes(x=method,y=diff)) +  geom_boxplot()+ ylim(-1,2)   

## W estimation evaluation
cbind(true_w[[2]][,1],true_w[[1]][,1],original_X[[1]]$w_tilde[,1],W_tilde[,1],w_est_alpha05[[1]][,1],w_est_alpha1[[1]][,1],w_est_alpha0[[1]][,1])

1*rbinom(1,1,0.4)+rnorm(1,0,1)*rbinom(1,1,0.8)

i=1
plot(w_est_alpha05[[i]],true_w[[i]], xlab = 'estmation',ylab='true', main="w_est_alpha05 t=1")
plot(w_est_alpha1[[i]],true_w[[i]], xlab = 'estmation',ylab='true', main="w_est_alpha1 t=1")
plot(original_X[[i]]$w_tilde,true_w[[i]], xlab = 'estmation',ylab='true', main="original w t=1")
plot(W_tilde,true_w[[i]], xlab = 'W tilde',ylab='true', main="W_tilde t=1")

compare.w = data.frame(diff = c(as.vector(w_est_alpha05[[i]]-true_w[[i]]),
                                as.vector(w_est_alpha1[[i]]-true_w[[i]]),
                                as.vector(w_est_alpha0[[i]]-true_w[[i]]), 
                                as.vector(W_tilde-true_w[[i]]), 
                                as.vector(original_X[[i]]$w_tilde-true_w[[i]])), 
                     method = c(rep("alpha_0.5",K*D),rep("alpha_1",K*D),rep("alpha_0",K*D),rep("W_tilde",K*D),rep("W_tilde_orig",K*D)))
#boxplot(diff~method, data=compare.w)
ggplot(compare.w, aes(x=method,y=diff)) +  geom_boxplot()+ ylim(-1,2) 

##


## evaluate gamma
hist(rst05$theta2$pi_ber)
hist(rst1$theta2$pi_ber)            

#z_est <-   rst$theta1$z 
#w_est <-  rst$theta1$w  
# tau_e_est <-  rst$theta2$tau_e  
# alpha_unif_est  <- rst$theta2$alpha_unif  
# gamma_est  <-  rst$theta2$gamma  
# v_est  <-   rst$theta2$v 
# pi_ber_est  <-  rst$theta2$pi_ber  
# tau_x_est  <-  rst$theta2$tau_x  
# tau_w_est  <-   rst$theta2$tau_w 
i=1
average_gamma_est_05 <- Reduce('+',rst05$theta2$gamma)/mcmc_samples_theta1
tbl_05 <- data.frame('Estmated_Gamma' = average_gamma_est_05[,i],"True_Gamma" = true_gamma[,i])
tbl_05$corruption_prob <- factor(corrupt_pi)
tbl_05$beta <- tau_xd_beta_para
tbl_05$Estimated_Gamma_bin <- 1*(tbl_05$Estmated_Gamma >= .5)
tbl_05$alpha <- 0.5

average_gamma_est_1 <- Reduce('+',rst1$theta2$gamma)/mcmc_samples_theta1
tbl_1 <- data.frame('Estmated_Gamma' = average_gamma_est_1[,i],"True_Gamma" = true_gamma[,i])
tbl_1$corruption_prob <- factor(corrupt_pi)
tbl_1$beta <- tau_xd_beta_para
tbl_1$Estimated_Gamma_bin <- 1*(tbl_1$Estmated_Gamma >= .5)
tbl_1$alpha <- 1

average_gamma_est_0 <- Reduce('+',rst0$theta2$gamma)/mcmc_samples_theta1
tbl_0 <- data.frame('Estmated_Gamma' = average_gamma_est_0[,i],"True_Gamma" = true_gamma[,i])
tbl_0$corruption_prob <- factor(corrupt_pi)
tbl_0$beta <- tau_xd_beta_para
tbl_0$Estimated_Gamma_bin <- 1*(tbl_0$Estmated_Gamma >= .5)
tbl_0$alpha <- 0

print(table(tbl_05$Estimated_Gamma_bin,tbl_05$True_Gamma))
print(table(tbl_1$Estimated_Gamma_bin,tbl_1$True_Gamma))

tbl <- rbind(tbl_05,tbl_1,tbl_0)

write.table(tbl,'./write/tbl.tauXdBeta=0.5.tauW.txt',quote = F, row.names = F)
## for debug
#hist((v_est[[k]]*gamma_est[[k]])[which(gamma_est[[k]]==1)],100)



# chisq.test(tbl)
# plot(v_est[[k]][which(true_gamma[,i]==1),i],true_v[which(true_gamma[,i]==1),i], xlab = 'estmation',ylab='true')
# plot(v_est[[k]][which(gamma_est[[k]][,i]==1),i],true_v[which(gamma_est[[k]][,i]==1),i], xlab = 'estmation',ylab='true', main='v')
# cor(v_est[[k]][which(gamma_est[[k]][,i]==1),i],true_v[which(gamma_est[[k]][,i]==1),i])



## ROC curve
ggplot(tbl, aes(d=True_Gamma,m=Estmated_Gamma,color=factor(alpha)))+geom_roc(n.cuts = 0)
ggplot(tbl, aes(d=True_Gamma,m=Estmated_Gamma,color=factor(tauW)))+geom_roc(n.cuts = 0)
ggplot(tbl, aes(d=True_Gamma,m=Estmated_Gamma,color=corruption_prob))+labs(fill="corruption_prob")+geom_roc(n.cuts = 0)





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
compare = data.frame(diff = c(as.vector(z_est_alpha0-true_z),as.vector(z_est_alpha1-true_z), as.vector(z_est_alpha05-true_z), as.vector(z_est_nnls-true_z),
                              as.vector(z_est_music-true_z)), 
                     method = c(rep("alpha_0",K*N),rep("alpha_1",K*N),rep("alpha_0.5",K*N), rep("NNLS",K*N),rep("MuSiC",K*N)),
                     para = "simulation11")                                          
#write.table(compare,'./write/simulation11.txt',quote = F,sep='\t',row.names = F)
ggplot(compare, aes(x=method,y=diff)) +  geom_boxplot()+ ylim(-.3,.4) 
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
