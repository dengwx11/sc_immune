set.seed(2021)
library(Seurat)
library(ggplot2)
library(reshape2)
library(Ckmeans.1d.dp)
source('Batch_Correction/run_liger.R')

## input demo
taget_tissue <- 'PBMC'
files = list.files(path = '/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk',pattern = "*_updated.rds", full.names = TRUE)
tissue_list <- c("BM","CB","Kidney","Liver","Lung","PBMC")
celltype_list <- c("B", "Neutrophil", "NK cell", "T")
YSG <- readRDS("/gpfs/ysm/pi/zhao-data/wd262/sc_immune/sc_immune/data/NSCLC/sg.list.rds")
output_path = "/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/pipeline_on_HCL"



get_tissue_specific_input <- function(taget_tissue,tissue_list,celltype_list,files,YSG,output_path){
    liger_output <- run_liger(files,tissue_list,YSG,output_path)
    seur.TPM_list <- liger_output$seur.TPM_list
    seur_list <- liger_output$seur_list
    liger <- liger_output$liger
    YSG <- intersect(YSG,liger@var.genes)
    W_tilde <- get_WH_W_tilde(liger,seur.TPM_list,seur_list)
    target_idx <- which(tissue_list == taget_tissue)
    tissue_gene_list <- get_tissue_gene(target_idx, W_tilde,tissue_list)
    w_tissue_indicator <- output_tissue_indicator(tissue_gene_list)
    w_tissue_indicator <- integrate_celltype(w_tissue_indicator,celltype_list,YSG)
    
    ans <- list()
    ans$YSG <- YSG
    ans$w_tissue_indicator <- w_tissue_indicator
    ans$X <- seur_list
    ans$seur.TPM_list <- seur.TPM_list
    ans$liger <- liger
    return(ans)
    
}