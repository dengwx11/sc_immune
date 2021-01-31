library(rliger)
library(Seurat)
library(umap)
library(ggplot2)
set.seed(2021)

## Data Loading and tramsform to log(TPM+1)
# NSCLC_seur <- readRDS("data/NSCLC/NSCLC_seur.rds") ## TPM
# NSCLC_seur[['RNA']]@data <- log1p(NSCLC_seur[['RNA']]@counts)
# NSCLC_seur.logTPM <- GetAssayData(NSCLC_seur[['RNA']],slot='data')

# PBMC_seur <- readRDS("data/NSCLC/PBMC_ImmuneCell_updated.rds") ## counts
# PBMC_seur <- NormalizeData(PBMC_seur,normalization.method = "LogNormalize",scale.factor = 1000000)
# PBMC_seur.logTPM <- GetAssayData(PBMC_seur[['RNA']],slot='data')


# Lung_seur <- readRDS("data/NSCLC/Lung_ImmuneCell_updated.rds") ## counts
# Lung_seur <- NormalizeData(Lung_seur,normalization.method = "LogNormalize",scale.factor = 1000000)
# Lung_seur.logTPM <- GetAssayData(Lung_seur[['RNA']],slot='data')

## Data Loading and tramsform to TPM
NSCLC_seur <- readRDS("data/NSCLC/NSCLC_seur.rds") ## TPM
NSCLC_seur[['RNA']]@data <- log1p(NSCLC_seur[['RNA']]@counts)
NSCLC_seur.TPM <- GetAssayData(NSCLC_seur[['RNA']],slot='counts')
NSCLC_seur <- NormalizeData(NSCLC_seur)
NSCLC_seur <- FindVariableFeatures(NSCLC_seur, selection.method = "vst", nfeatures = 500)
NSCLC_seur <- ScaleData(NSCLC_seur)
NSCLC_seur<- RunPCA(NSCLC_seur, verbose = FALSE)
NSCLC_seur <- FindNeighbors(NSCLC_seur, dims = 1:25)
NSCLC_seur <- FindClusters(NSCLC_seur, resolution = 0.8)
NSCLC_seur <- RunUMAP(NSCLC_seur, dims = 1:50)

PBMC_seur <- readRDS("data/NSCLC/PBMC_ImmuneCell_updated.rds") ## counts
PBMC_seur <- NormalizeData(PBMC_seur,normalization.method = "LogNormalize",scale.factor = 1000000)
PBMC_seur.TPM <- exp(GetAssayData(PBMC_seur[['RNA']],slot='data'))-1
PBMC_seur <- NormalizeData(PBMC_seur)
PBMC_seur <- FindVariableFeatures(PBMC_seur, selection.method = "vst", nfeatures = 500)
PBMC_seur <- ScaleData(PBMC_seur)
PBMC_seur<- RunPCA(PBMC_seur, verbose = FALSE)
PBMC_seur <- FindNeighbors(PBMC_seur, dims = 1:25)
PBMC_seur <- FindClusters(PBMC_seur, resolution = 0.8)
PBMC_seur <- RunUMAP(PBMC_seur, dims = 1:50)


Lung_seur <- readRDS("data/NSCLC/Lung_ImmuneCell_updated.rds") ## counts
Lung_seur <- NormalizeData(Lung_seur,normalization.method = "LogNormalize",scale.factor = 1000000)
Lung_seur.TPM <- exp(GetAssayData(Lung_seur[['RNA']],slot='data'))-1
Lung_seur <- NormalizeData(Lung_seur)
Lung_seur <- FindVariableFeatures(Lung_seur, selection.method = "vst", nfeatures = 500)
Lung_seur <- ScaleData(Lung_seur)
Lung_seur<- RunPCA(Lung_seur, verbose = FALSE)
Lung_seur <- FindNeighbors(Lung_seur, dims = 1:25)
Lung_seur <- FindClusters(Lung_seur, resolution = 0.8)
Lung_seur <- RunUMAP(Lung_seur, dims = 1:50)

YSG <- readRDS("data/NSCLC/sg.list.rds")

## Stage I: preprocessing and normalization
# ifnb_liger <- createLiger(list(NSCLC = NSCLC_seur.logTPM, 
#                                PBMC_HCL = PBMC_seur.logTPM))
# ifnb_liger <- createLiger(list(NSCLC = NSCLC_seur.TPM/100, 
#                                PBMC_HCL = PBMC_seur.TPM/100))
ifnb_liger <- createLiger(list(NSCLC = NSCLC_seur.TPM/100, 
                               PBMC_HCL = PBMC_seur.TPM/100,
                               Lung_HCL = Lung_seur.TPM/100))                               
ifnb_liger <- normalize(ifnb_liger)
ifnb_liger <- selectGenes(ifnb_liger)
ifnb_liger@var.genes <- intersect(YSG,ifnb_liger@var.genes)
ifnb_liger <- scaleNotCenter(ifnb_liger) 

## Stage II: joint matrix factorization
ifnb_liger <- optimizeALS(ifnb_liger, k = 20)
saveRDS(ifnb_liger,"write/liger/NSCLC_PBMC_Lung_TPM.rds")

## Stage III: quantile normalization and joint clustering
ifnb_liger <- quantile_norm(ifnb_liger)
ifnb_liger <- louvainCluster(ifnb_liger, resolution = 0.25)

## Stage IV: visualization and downstream analysis
ifnb_liger <- runTSNE(ifnb_liger) 
ifnb_liger <- runUMAP(ifnb_liger)
saveRDS(ifnb_liger,"write/liger/NSCLC_PBMC_Lung_TPM.rds")
ifnb_liger <- runTSNE(ifnb_liger,use.raw=TRUE) 
ifnb_liger <- runUMAP(ifnb_liger,use.raw=TRUE)
saveRDS(ifnb_liger,"write/liger/NSCLC_PBMC_Lung_TPM_raw.rds")



## Stage IV: visualization and downstream analysis
ifnb_liger <- readRDS("write/liger/NSCLC_PBMC_Lung_TPM.rds")
p <- plotByDatasetAndCluster(ifnb_liger, return.plots = T)
plot_grid(p[[1]], p[[2]])

plotByDatasetAndCluster(ifnb_liger, axis.labels = c('UMAP 1', 'UMAP 2'))

gene_loadings <- plotGeneLoadings(ifnb_liger, do.spec.plot = FALSE, return.plots = TRUE)
gene_loadings[[13]]

plotGene(ifnb_liger, "CD3E")
plotGene(ifnb_liger, "CD4")
plotGene(ifnb_liger, "CD8A")
plotGene(ifnb_liger, "CD14")
plotGene(ifnb_liger, "CD19")
plotGene(ifnb_liger, "PTPRC")

plotClusterProportions(ifnb_liger)
plotClusterFactors(ifnb_liger, use.aligned = T)

mydir = "plots/liger/NSCLC_PBMC_Lung_"
pdf(paste0(mydir, "word_clouds.pdf"))
plotWordClouds(ifnb_liger, dataset1 = "NSCLC", dataset2 = "PBMC_HCL")
dev.off()

pdf(paste0(mydir, "plot_factors.pdf"))
plotFactors(ifnb_liger)
dev.off()

markers <- getFactorMarkers(ifnb_liger, dataset1 = "NSCLC", dataset2 = "PBMC_HCL", num.genes = 10)
plotGene(ifnb_liger, gene = "CD14")
plotGeneViolin(ifnb_liger, gene = "CD14")



## batch correction
WH_NSCLC <- t(ifnb_liger@W)%*%t(ifnb_liger@H.norm[colnames(NSCLC_seur.TPM),])
#WH_NSCLC <- t(ifnb_liger@W)%*%t(ifnb_liger@H$NSCLC)
WVH_NSCLC <- (t(ifnb_liger@W)+t(ifnb_liger@V$NSCLC))%*%t(ifnb_liger@H.norm[colnames(NSCLC_seur.TPM),])
#WVH_NSCLC <- (t(ifnb_liger@W)+t(ifnb_liger@V$NSCLC))%*%t(ifnb_liger@H$NSCLC)
mat_NSCLC <- GetAssayData(NSCLC_seur,slot="counts")[ifnb_liger@var.genes,]
mat_NSCLC <- as.matrix(mat_NSCLC)
idx.sample <- sample(seq(dim(mat_NSCLC)[1]*dim(mat_NSCLC)[2]),1000)
plot(as.vector(WVH_NSCLC)[idx.sample],as.vector(mat_NSCLC)[idx.sample]/100)
plot(as.vector(WH_NSCLC)[idx.sample],as.vector(mat_NSCLC)[idx.sample]/100)

cor(as.vector(WVH_NSCLC)[idx.sample],as.vector(mat_NSCLC)[idx.sample]/100)
cor(as.vector(WH_NSCLC)[idx.sample],as.vector(mat_NSCLC)[idx.sample]/100)

WH_PBMC <- t(ifnb_liger@W)%*%t(ifnb_liger@H.norm[colnames(PBMC_seur.TPM),])
#WH_PBMC <- t(ifnb_liger@W)%*%t(ifnb_liger@H$PBMC_HCL)
WVH_PBMC<- (t(ifnb_liger@W)+t(ifnb_liger@V$PBMC_HCL))%*%t(ifnb_liger@H.norm[colnames(PBMC_seur.TPM),])
#WVH_PBMC<- (t(ifnb_liger@W)+t(ifnb_liger@V$PBMC_HCL))%*%t(ifnb_liger@H$PBMC_HCL)
mat_PBMC <- PBMC_seur.TPM[ifnb_liger@var.genes,]
mat_PBMC <- as.matrix(mat_PBMC)
idx.sample <- sample(seq(dim(mat_PBMC)[1]*dim(mat_PBMC)[2]),2000)
plot(as.vector(WVH_PBMC)[idx.sample],as.vector(mat_PBMC)[idx.sample])
plot(as.vector(WH_PBMC)[idx.sample],as.vector(mat_PBMC)[idx.sample])

cor(as.vector(WVH_PBMC)[idx.sample],as.vector(mat_PBMC)[idx.sample])
cor(as.vector(WH_PBMC)[idx.sample],as.vector(mat_PBMC)[idx.sample])

WH_Lung <- t(ifnb_liger@W)%*%t(ifnb_liger@H.norm[colnames(Lung_seur.TPM),])
#WH_PBMC <- t(ifnb_liger@W)%*%t(ifnb_liger@H$PBMC_HCL)
WVH_Lung<- (t(ifnb_liger@W)+t(ifnb_liger@V$Lung_HCL))%*%t(ifnb_liger@H.norm[colnames(Lung_seur.TPM),])
#WVH_PBMC<- (t(ifnb_liger@W)+t(ifnb_liger@V$PBMC_HCL))%*%t(ifnb_liger@H$PBMC_HCL)
mat_Lung <- Lung_seur.TPM[ifnb_liger@var.genes,]
mat_Lung <- as.matrix(mat_Lung)
idx.sample <- sample(seq(dim(mat_Lung)[1]*dim(mat_Lung)[2]),2000)
plot(as.vector(WVH_Lung)[idx.sample],as.vector(mat_Lung)[idx.sample])
plot(as.vector(WH_Lung)[idx.sample],as.vector(mat_Lung)[idx.sample])

cor(as.vector(WVH_Lung)[idx.sample],as.vector(mat_Lung)[idx.sample])
cor(as.vector(WH_Lung)[idx.sample],as.vector(mat_Lung)[idx.sample])

## How are zero expressions corrected on WH
idx <- which(as.vector(mat_NSCLC)<0.01)
idx.sample <- sample(seq(length(idx)),1000)
hist(as.vector(WH_NSCLC)[idx[idx.sample]],100,main="Histogram of WH (corrected expression) when the raw counts<0.01")

idx <- which(as.vector(mat_PBMC)<0.01)
idx.sample <- sample(seq(length(idx)),1000)
hist(as.vector(WH_PBMC)[idx[idx.sample]],100,main="Histogram of WH (corrected expression) when the raw counts<0.01")

## visualization on W*H.norm
new_mat <- cbind(WH_NSCLC,WH_PBMC,WH_Lung)
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
saveRDS(new_seur,'write/liger/BatchCorrected_NSCLC_PBMC_Lung_seur.rds')


new_split <- SplitObject(new_seur,split.by = 'Batch')
new_NSCLC <- new_split[['NSCLC']]
new_PBMC <- new_split[['PBMC_HCL']]
new_Lung <- new_split[['Lung_HCL']]

p1 <- DimPlot( new_NSCLC , group.by='Celltype_used')+ ggtitle("NSCLC, Batch Corrected")
p2 <- DimPlot( NSCLC_seur , group.by='Celltype_used')+ ggtitle("NSCLC, Raw")
plot_grid(p1, p2)

p1 <- DimPlot( new_PBMC , group.by='Celltype_used') + ggtitle("PBMC HCL, Batch Corrected")
p2 <- DimPlot( PBMC_seur , group.by='Celltype_used') + ggtitle("PBMC HCL, Raw")
plot_grid(p1, p2)

p1 <- DimPlot( new_Lung , group.by='Celltype_used') + ggtitle("Lung HCL, Batch Corrected")
p2 <- DimPlot( Lung_seur , group.by='Celltype_used') + ggtitle("Lung HCL, Raw")
plot_grid(p1, p2)

## save corrected and normalized data in Seurat Object data slot
NSCLC_seur@assays$RNA@data <- WH_NSCLC
PBMC_seur@assays$RNA@data <- WH_PBMC
Lung_seur@assays$RNA@data <- WH_Lung

celltype.list <- c("B", "Monocyte", "NK cell", "T")
new_PBMC <- subset(new_PBMC, Celltype_used %in% celltype.list)
set.seed(2021)
idx.sample <- sample(seq(ncol(new_PBMC)),4000)
new_PBMC <- subset(new_PBMC, cells=Cells(new_PBMC)[idx.sample])
new_Lung <- subset(new_Lung, Celltype_used %in% celltype.list)
saveRDS(new_NSCLC,"data/NSCLC/NSCLC_liger_seur.rds")
saveRDS(new_PBMC,"data/NSCLC/PBMC_liger_seur.rds")
saveRDS(new_Lung,"data/NSCLC/Lung_liger_seur.rds")

