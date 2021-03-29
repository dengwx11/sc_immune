options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

T <- as.numeric(args[1])
D <- as.numeric(args[2])
K <- as.numeric(args[3])
corrupt_pi <- as.numeric(args[4])
tau_w_para <- as.numeric(args[5])
tau_xd_beta_para <- as.numeric(args[6])
mcmc_samples_theta1 <- as.numeric(args[7])
seed <- as.integer(args[8])

set.seed(seed)
library(nnls)
library(MuSiC)
library(xbioc)
library(ggplot2)
library(Biobase)
library(plotROC)
#set.seed(2021)
source('SAME/run_same.v2.R')
source('Batch_Correction/platform_batch_correction_combat.R')

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
Lambda = c(0:mcmc_samples_theta1) # Lambda = c(0,1,2,3,...,100)

target_tissue <- 'PBMC'
files = list.files(path = '/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk',pattern = "*_updated.rds", full.names = TRUE)
tissue_list <- c("BM","CB","Kidney","Liver","Lung","PBMC")
celltype_used_list <- c("B", "Neutrophil", "NK cell", "T")
YSG <- readRDS("/gpfs/ysm/pi/zhao-data/wd262/sc_immune/sc_immune/data/NSCLC/sg.list.rds")
output_path = "/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/pipeline_on_HCL"
mcmc_samples_theta1=5
Y0 <- read.table("data/NSCLC/Fig2b-WholeBlood_RNAseq.txt", row.names = 1, header = T, sep = "\t")

rst <- run_SAME(target_tissue,tissue_list,celltype_used_list,files,YSG,empirical_pi=0.3,mcmc_samples_theta1,
                liger.turnon=FALSE,output_path)
# rst_bulk <- Y_batch_correct(Y0, rst$X[[target_tissue]], rst$YSG)
rst_frac <- fraction_est(Y0, rst,target_tissue, celltype_used_list)

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


rst <- list()
rst$same_rst <- rst1
rst$same_input <- same_input
rst$music <- z_est_music
rst$nnls_music <- z_est_nnls
mydir <- "/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/simulation_rst/"            
saveRDS(rst,paste0(mydir,"rst.alpha=",alpha,'.',str_para,'.rds'))