library(Seurat)

#read dataset
seur <- readRDS("/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/Lung_ImmuneCell_updated.rds")

#output file
PseudoBulk_Lung <- data.frame()

mat <- seur@assays$RNA@counts
meta_df <- seur[["Celltype_used"]]
Lung_stat <- as.data.frame(table(meta_df))
Lung_stat[,2] = Lung_stat[,2]/sum(Lung_stat[,2])
set.seed(2020)
times = 50
for(i in 1:times){
	cells <- sample(Cells(seur), 10000, replace = T)
	mat_new <- mat[,cells]
	meta_df_new <- meta_df[cells,]
	stat <- as.data.frame(table(meta_df_new))
	Lung_stat[,i+2] =stat[,2]/10000
	for (j in 1:length(rownames(mat_new))) {
		PseudoBulk_Lung[j,i] = round(sum(mat_new[j,]) + var(mat_new[j,])/50)
	}
}
rownames(Lung_stat) = Lung_stat[,1]
Lung_stat = Lung_stat[,-1]
saveRDS(PseudoBulk_Lung, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/PseudoBulk_Lung.rds")
saveRDS(Lung_stat, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/Lung_stat.rds")


