source('SAME/run_same.R')
ifnb_liger <- readRDS("write/liger/NSCLC_PBMC_Lung_TPM.rds")


WVH_NSCLC <- (t(ifnb_liger@W)+t(ifnb_liger@V$PBMC_HCL))%*%t(ifnb_liger@H.norm[colnames(NSCLC_seur.TPM),])

WVH_PBMC<- (t(ifnb_liger@W)+t(ifnb_liger@V$PBMC_HCL))%*%t(ifnb_liger@H.norm[colnames(PBMC_seur.TPM),])

WVH_Lung<- (t(ifnb_liger@W)+t(ifnb_liger@V$PBMC_HCL))%*%t(ifnb_liger@H.norm[colnames(Lung_seur.TPM),])


## visualization on W*H.norm
new_mat <- cbind(WVH_NSCLC,WVH_PBMC,WVH_Lung)
Celltype_used <-  c(as.character(NSCLC_seur$Celltype_used), 
            as.character(PBMC_seur$Celltype_used), 
            as.character(Lung_seur$Celltype_used))
Batch <- c(rep("NSCLC",ncol(NSCLC_seur)),rep("PBMC_HCL",ncol(PBMC_seur)),rep("Lung_HCL",ncol(Lung_seur)))
new_seur <- CreateSeuratObject(new_mat)
new_seur$Celltype_used <- Celltype_used
new_seur$Batch <- Batch

## umap visualization
# standard log-normalization
new_seur <- NormalizeData(new_seur)
new_seur@assays$RNA@data <- new_seur@assays$RNA@counts
#new_seur@assays$RNA@data <- log1p(new_seur@assays$RNA@counts)
# choose ~1k variable features
new_seur <- FindVariableFeatures(new_seur, selection.method = "vst", nfeatures = 500)
# standard scaling (no regression)
new_seur <- ScaleData(new_seur)
# Run PCA, select 13 PCs for tSNE visualization and graph-based clustering
new_seur<- RunPCA(new_seur, verbose = FALSE)

new_seur <- FindNeighbors(new_seur, dims = 1:25)
new_seur <- FindClusters(new_seur, resolution = 0.8)
new_seur <- RunUMAP(new_seur, dims = 1:50)
DimPlot(new_seur,group.by = 'Batch')
DimPlot(new_seur,group.by = 'Celltype_used')

Y0 <- read.table("data/NSCLC/Fig2b-WholeBlood_RNAseq.txt", row.names = 1, header = T, sep = "\t")
Y0_sc <- as.matrix(Y0[rownames(new_seur),])
X <- list()

new_split <- SplitObject(new_seur,split.by = 'Batch')
new_NSCLC <- new_split[['NSCLC']]
new_PBMC <- new_split[['PBMC_HCL']]
new_Lung <- new_split[['Lung_HCL']]

X <- list()
X[[1]] <- subset(new_PBMC, Celltype_used %in% celltype.list)
X[[2]] <- new_NSCLC
X[[3]] <- subset(new_Lung, Celltype_used %in% celltype.list)


Idents(X[[1]]) <- "Celltype_used"
X[[1]] <- NormalizeData(X[[1]])
W_tilde <- as.matrix(AverageExpression(X[[1]],slot='data')$RNA)/100
Idents(X[[2]]) <- "Celltype_used"
X[[2]] <- NormalizeData(X[[2]])
W_tilde_2 <- as.matrix(AverageExpression(X[[2]],slot='data')$RNA)/100
Idents(X[[3]]) <- "Celltype_used"
X[[3]] <- NormalizeData(X[[3]])
W_tilde_3 <- as.matrix(AverageExpression(X[[3]],slot='data')$RNA)/100
## SAME
mcmc_samples_theta1 = 30
Lambda = c(0:mcmc_samples_theta1)
c_k <- lapply(c(1:2), function(i) table(X[[i]]$Celltype_used))
c_k <- t(Reduce(rbind,c_k))
YSG <- rownames(new_seur)

T=length(X)
D=nrow(Y0_sc)
K=ncol(W_tilde)
N = ncol(Y0)

rst1 <- SAME(Y0_sc, X, W_tilde,
            mcmc_samples_theta1, Lambda, c_k, YSG, alpha =1)
w_est_alpha1 <-  rst1$theta1$w
est_vg <- lapply(c(1:mcmc_samples_theta1), function(i) rst1$theta2$gamma[[i]] * rst1$theta2$v[[i]])
average_est_vg_1 <- Reduce('+',est_vg)/mcmc_samples_theta1 

celltype.list <- c("B", "Monocyte", "NK cell", "T")
cor(average_est_vg_1,W_tilde[,celltype.list])
plot(average_est_vg_1,W_tilde[,celltype.list])

celltype.intersect <- intersect(colnames(W_tilde_2),colnames(W_tilde))
plot(W_tilde_2[,celltype.intersect],W_tilde_3[,celltype.intersect])
cor(W_tilde[,celltype.intersect],W_tilde_2[,celltype.intersect])

library(reshape2)
W_tilde_A <- W_tilde
W_tilde_B <- W_tilde_3
celltype.intersect <- intersect(colnames(W_tilde_A),colnames(W_tilde_B))
df.cor1 <- melt(W_tilde_A[,celltype.intersect],value.name = "PBMC_HCL")
df.cor2 <- melt(W_tilde_B[,celltype.intersect],value.name = "Lung_HCL")
df.cor <- cbind(df.cor1,Lung_HCL=df.cor2[,3])
ggplot(df.cor,aes(x=PBMC_HCL,y=Lung_HCL))+geom_point(aes(colour = factor(Var2)),size=.1)