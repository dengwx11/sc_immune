library(Seurat)

#read dataset
seur <- readRDS("/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/Liver_ImmuneCell_updated.rds")

#output file
PseudoBulk_Liver <- data.frame()

mat <- seur@assays$RNA@counts
meta_df <- seur[["Celltype_used"]]
Liver_stat <- as.data.frame(table(meta_df))
Liver_stat[,2] = Liver_stat[,2]/sum(Liver_stat[,2])
set.seed(2020)
times = 50
for(i in 1:times){
	cells <- sample(Cells(seur), 10000, replace = T)
	mat_new <- mat[,cells]
	meta_df_new <- meta_df[cells,]
	stat <- as.data.frame(table(meta_df_new))
	Liver_stat[,i+2] =stat[,2]/10000
	for (j in 1:length(rownames(mat_new))) {
		PseudoBulk_Liver[j,i] = round(sum(mat_new[j,]) + var(mat_new[j,])/50)
	}
}
rownames(Liver_stat) = Liver_stat[,1]
Liver_stat = Liver_stat[,-1]
saveRDS(PseudoBulk_Liver, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/PseudoBulk_Liver.rds")
saveRDS(Liver_stat, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/Liver_stat.rds")


