library(Seurat)

#read dataset
seur <- readRDS("/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/BM_ImmuneCell_updated.rds")

#output file
PseudoBulk_BM <- data.frame()
BM_stat <- data.frame()

set.seed(2020)
times = 50
for(i in 1:times){
	cells <- sample(Cells(seur), 10000, replace = T)
	mat <- seur@assays$RNA@counts
	mat_new <- mat[,cells]
	meta_df <- seur[["Celltype_used"]]
	meta_df_new <- meta_df[cells,]
	stat <- as.data.frame(table(meta_df_new))
	BM_stat = stat
	BM_stat[,i+1] =stat[,2]
	for (j in 1:length(rownames(mat_new))) {
		PseudoBulk_Lung[j,i] = round(sum(mat_new[j,]) + var(mat_new[j,])/50)
	}
}
rownames(BM_stat) = BM_stat[,1]
BM_stat = BM_stat[,-1]/10000
saveRDS(PseudoBulk_BM, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/PseudoBulk_BM.rds")
saveRDS(BM_stat, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/BM_stat.rds")

