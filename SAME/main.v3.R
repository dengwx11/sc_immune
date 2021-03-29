setwd("~/zhao-data/sc_immune/sc_immune/")
source('SAME/run_same.v2.R')
source('Batch_Correction/platform_batch_correction_combat.R')

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
rst_frac <- fraction_est(Y0, rst,target_tissue, celltype_used_list,method = c("tranSig","music","w_empirical"), adj.to.sc = T)