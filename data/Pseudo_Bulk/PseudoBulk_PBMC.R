library(Seurat)

#read dataset

seur <- readRDS("/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/PBMC_ImmuneCell_updated.rds")


#output file
PseudoBulk_PBMC <- data.frame()

mat <- seur@assays$RNA@counts
meta_df <- seur[["Celltype_used"]]
PBMC_stat <- as.data.frame(table(meta_df))
PBMC_stat[,2] = PBMC_stat[,2]/sum(PBMC_stat[,2])
set.seed(2020)
times = 50
for(i in 1:times){
	cells <- sample(Cells(seur), 10000, replace = T)
	mat_new <- mat[,cells]
	meta_df_new <- meta_df[cells,]
	stat <- as.data.frame(table(meta_df_new))
	PBMC_stat[,i+2] =stat[,2]/10000
	for (j in 1:length(rownames(mat_new))) {
		PseudoBulk_PBMC[j,i] = round(sum(mat_new[j,]) + var(mat_new[j,])/50)
	}
}
rownames(PBMC_stat) = PBMC_stat[,1]
PBMC_stat = PBMC_stat[,-1]
saveRDS(PseudoBulk_PBMC, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/PseudoBulk_PBMC.rds")
saveRDS(PBMC_stat, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/PBMC_stat.rds")

