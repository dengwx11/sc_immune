set.seed(2021)
library(Seurat)
library(ggplot2)
library(reshape2)
library(Ckmeans.1d.dp)
library(rliger)
setwd("/data4/lbl/Yale")
source("SAME/run_same.v2.R")
source('Batch_Correction/run_liger.R')
source('SAME/platform_batch_correction.R')

## input demo
taget_tissue <- 'PBMC'
# files = list.files(path = '/data4/lbl/Yale/HCL/protein_coding',pattern = "*_updated.rds", full.names = TRUE)
tissue_list <- c("PBMC","Lung","BM","CB"
                 ,"Kidney","Liver"
)
files = sapply(1:6, function(i)paste0('/data4/lbl/Yale/HCL/protein_coding/',tissue_list[i], "_ImmuneCell_updated.rds"))
celltype.list <- c("B", "Neutrophil", "NK cell", "T")

YSG <- readRDS("/data4/lbl/Yale/HCL/SAME/sg.list (1).rds")
output_path = "/data4/lbl/Yale/SAME/Batch_Correction/"

#Get w_tissue_indicator
ans <- get_tissue_specific_input(taget_tissue,tissue_list,celltype.list,files,YSG,output_path)
w_tissue_indicator <- ans$w_tissue_indicator
YSG <- ans$YSG
X <- ans$X

# clean format
table(X[[1]]$Celltype_used)
for (i in 1:length(tissue_list)) {
  seur <- X[[i]]
  X[[i]] <- subset(seur, Celltype_used %in% celltype.list)
}
table(X[[1]]$Celltype_used)
names(w_tissue_indicator)
Cl <- get_Cl(X)
print(Cl$levels)

#prepare SAME input
T = length(X)
K = length(Cl$levels)
c_k = matrix(0, nrow = K, ncol = T)
D = length(YSG)
for (t in 1:T){
  seur = X[[t]]
  for (k in 1:K){
    seur.list <- SplitObject(seur, split.by = "Celltype_used")
    if (is.null(seur.list[[celltype.list[k]]])){
      c_k[k,t] = 0
    } else {
      c_k[k,t] = ncol(seur.list[[celltype.list[k]]])
    }
  }
}
C0 = sum(c_k)

mcmc_samples_theta1 = 50
Lambda = c(0:mcmc_samples_theta1)


#read bulk seq and adjust into single cell space
Y0 <- read.table("/data4/lbl/Yale/HCL/Bayesian/CibersortX/Fig2ab-NSCLC_PBMCs/Fig2b-WholeBlood_RNAseq.txt", row.names = 1, header = T, sep = "\t")
YSG <- intersect(YSG, rownames(Y0))
rst <- Y_batch_correct(Y0, X[[1]], YSG)
Y0_sc <- rst$bulk_to_sc
Y0_sc <- Y0_sc[YSG,] #TPM space
N = ncol(Y0_sc)

for(i in 1:T){
  w_indicator <- w_tissue_indicator[[i]]
  w_tissue_indicator[[i]] <- w_indicator[YSG,]
}

#read ground truth
WBC_gt <- read.table("/data4/lbl/Yale/HCL/Bayesian/CibersortX/Fig2ab-NSCLC_PBMCs/Fig2b_ground_truth_whole_blood.txt", row.names = 1, header = 1, sep = "\t")
gt_sub <- WBC_gt[,c(7,1,8,4)]

# run SAME with a series of pi
rst_new_step0.1.list <- list()
corr_celltype <- c()
empirical_pi <- c(0.1, 0.2, 0.3, 0.4, 0.5,  0.6, 0.7,  0.8, 0.9)
for (i in 1:length(empirical_pi)) {
  rst_bacth_corrected <- SAME(X,w_tissue_indicator,empirical_pi[i],mcmc_samples_theta1, Lambda, c_k, YSG) #time  554.109s
  average_est_vg_1 <- rst_bacth_corrected$vg
  #estimate z by nnls
  z_est_vg_nnls <- t(sapply(c(1:N), function(j) nnls(average_est_vg_1,(Y0_sc[YSG,j]))$x))
  corr_vg_nnls<- sapply(c(1:4), function(i) cor(as.vector(as.matrix(gt_sub)[,i]), as.vector(z_est_vg_nnls[,i])))

  rst_bacth_corrected$theta2$corr_celltype <- corr_vg_nnls
  rst_new_step0.1.list[[i]] <- rst_bacth_corrected
  saveRDS(rst_bacth_corrected, file = paste0(output_path, "rst_new_tauv0.6_pi", empirical_pi[i], ".rds"))
}

save(rst_new_step0.1.list, file = "SAME/Batch_Correction/rst/rst_new_step0.1.list.RData")


