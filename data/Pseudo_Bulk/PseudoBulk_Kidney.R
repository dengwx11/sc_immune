library(Seurat)

#read dataset
seur <- readRDS("/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/Kidney_ImmuneCell_updated.rds")



#output file
PseudoBulk_Kidney <- data.frame()

mat <- seur@assays$RNA@counts
meta_df <- seur[["Celltype_used"]]
Kidney_stat <- as.data.frame(table(meta_df))
Kidney_stat[,2] = Kidney_stat[,2]/sum(Kidney_stat[,2])
set.seed(2020)
times = 50
for(i in 1:times){
	cells <- sample(Cells(seur), 10000, replace = T)
	mat_new <- mat[,cells]
	meta_df_new <- meta_df[cells,]
	stat <- as.data.frame(table(meta_df_new))
	Kidney_stat[,i+2] =stat[,2]/10000
	for (j in 1:length(rownames(mat_new))) {
		PseudoBulk_Kidney[j,i] = round(sum(mat_new[j,]) + var(mat_new[j,])/50)
	}
}
rownames(Kidney_stat) = Kidney_stat[,1]
Kidney_stat = Kidney_stat[,-1]
saveRDS(PseudoBulk_Kidney, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/PseudoBulk_Kidney.rds")
saveRDS(Kidney_stat, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/Kidney_stat.rds")


