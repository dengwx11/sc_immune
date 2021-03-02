set.seed(2021)

taget_tissue <- 'PBMC'


library(Seurat)
library(ggplot2)
library(reshape2)
library(Ckmeans.1d.dp)
source('Batch_Correction/run_liger.R')


files = list.files(path = '/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk',pattern = "*_updated.rds", full.names = TRUE)
tissue_list <- c("BM","CB","Kidney","Liver","Lung","PBMC")
YSG <- readRDS("/gpfs/ysm/pi/zhao-data/wd262/sc_immune/sc_immune/data/NSCLC/sg.list.rds")
output_path = "/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/pipeline_on_HCL"

liger_output <- run_liger(files,tissue_list,YSG,output_path)
seur.TPM_list <- liger_output$seur.TPM_list
seur_list <- liger_output$seur_list
liger <- liger_output$liger
YSG <- intersect(YSG,liger@var.genes)

W_tilde <- get_WH_W_tilde(liger,seur.TPM_list,seur_list)

target_idx <- which(tissue_list == taget_tissue)

tissue_gene_list <- get_tissue_gene(target_idx, W_tilde,tissue_list)
w_tissue_indicator <- output_tissue_indicator(tissue_gene_list)