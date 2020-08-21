library(Seurat)

#read dataset

seur <- readRDS("/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/CB_ImmuneCell_updated.rds")


#output file
PseudoBulk_CB <- data.frame()

mat <- seur@assays$RNA@counts
meta_df <- seur[["Celltype_used"]]
CB_stat <- as.data.frame(table(meta_df))
CB_stat[,2] = CB_stat[,2]/sum(CB_stat[,2])
set.seed(2020)
times = 50
for(i in 1:times){
	cells <- sample(Cells(seur), 10000, replace = T)
	mat_new <- mat[,cells]
	meta_df_new <- meta_df[cells,]
	stat <- as.data.frame(table(meta_df_new))
	CB_stat[,i+2] =stat[,2]/10000
	for (j in 1:length(rownames(mat_new))) {
		PseudoBulk_CB[j,i] = round(sum(mat_new[j,]) + var(mat_new[j,])/50)
	}
}
rownames(CB_stat) = CB_stat[,1]
CB_stat = CB_stat[,-1]
saveRDS(PseudoBulk_CB, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/PseudoBulk_CB.rds")
saveRDS(CB_stat, "/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk/CB_stat.rds")


