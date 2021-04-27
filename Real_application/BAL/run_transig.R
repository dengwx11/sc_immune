setwd("~/zhao-data/sc_immune/sc_immune/")
source('SAME/run_same.v2.R')
source('Batch_Correction/platform_batch_correction_combat.R')

celltype_used_lists <- readRDS('/gpfs/loomis/pi/zhao2/wd262/shared_data/bulk/E-MTAB-62/celltype_used_lists.rds')
file_tissue_list <-read.table('/gpfs/loomis/pi/zhao2/wd262/shared_data/bulk/E-MTAB-62/file_tissue_list.txt',header=T)
rownames(file_tissue_list) <- c(1:nrow(file_tissue_list))

target_tissue <- 'Adult-Lung'
idx <- c(1,3,5,9,14:16,19,26,29)
files = file_tissue_list$used_files[idx]
tissue_list <- file_tissue_list$used_tissue[idx]
celltype_used_list <- celltype_used_lists[['Adult-Lung']]
celltype_used_list <- celltype_used_list[which(celltype_used_list!='DELETED')]
YSG <- readRDS("/gpfs/ysm/pi/zhao-data/wd262/sc_immune/sc_immune/data/NSCLC/sg.list.rds")
output_path = "/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/real_application/BAL"
mcmc_samples_theta1=5
Y0 <- readRDS('/gpfs/ysm/pi/zhao-data/wd262/sc_immune/sc_immune/data/BAL/bulk.rawcount.rds')

rst <- run_SAME(target_tissue,tissue_list,celltype_used_list,files,YSG,empirical_pi=0.3,mcmc_samples_theta1,
                liger.turnon=TRUE,output_path)
saveRDS(rst,paste0(output_path,'/rst.rds'))
# rst_bulk <- Y_batch_correct(Y0, rst$X[[target_tissue]], rst$YSG)
rst_frac <- fraction_est(Y0, rst,target_tissue, celltype_used_list,method = c("tranSig","music","w_empirical"), adj.to.sc = T)
saveRDS(rst_frac,paste0(output_path,'/rst_frac.rds'))
