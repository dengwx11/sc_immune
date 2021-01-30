library(rliger)
library(Seurat)
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

PBMC_seur <- readRDS("data/NSCLC/PBMC_ImmuneCell_updated.rds") ## counts
PBMC_seur <- NormalizeData(PBMC_seur,normalization.method = "LogNormalize",scale.factor = 1000000)
PBMC_seur.TPM <- exp(GetAssayData(PBMC_seur[['RNA']],slot='data'))-1


Lung_seur <- readRDS("data/NSCLC/Lung_ImmuneCell_updated.rds") ## counts
Lung_seur <- NormalizeData(Lung_seur,normalization.method = "LogNormalize",scale.factor = 1000000)
Lung_seur.TPM <- exp(GetAssayData(Lung_seur[['RNA']],slot='data'))-1

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
ifnb_liger <- readRDS("write/liger/NSCLC_PBMC.HCL_TPM.rds")
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

mydir = "plots/liger/NSCLC_PBMC.HCL_"
pdf(paste0(mydir, "word_clouds.pdf"))
plotWordClouds(ifnb_liger, dataset1 = "NSCLC", dataset2 = "PBMC_HCL")
dev.off()

pdf(paste0(mydir, "plot_factors.pdf"))
plotFactors(ifnb_liger)
dev.off()

markers <- getFactorMarkers(ifnb_liger, dataset1 = "NSCLC", dataset2 = "PBMC_HCL", num.genes = 10)
plotGene(ifnb_liger, gene = "CD14")
plotGeneViolin(ifnb_liger, gene = "CD14")

cbind(ifnb_liger@H.norm[1:10,1],ifnb_liger@H$NSCLC[1:10,1])
sum(ifnb_liger@H.norm[3,])

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

## How are zero expressions corrected on WH
idx <- which(as.vector(mat_NSCLC)<0.01)
idx.sample <- sample(seq(length(idx)),1000)
hist(as.vector(WH_NSCLC)[idx[idx.sample]],100,main="Histogram of WH (corrected expression) when the raw counts<0.01")

idx <- which(as.vector(mat_PBMC)<0.01)
idx.sample <- sample(seq(length(idx)),1000)
hist(as.vector(WH_PBMC)[idx[idx.sample]],100,main="Histogram of WH (corrected expression) when the raw counts<0.01")

str(WH_NSCLC)
str(WH_PBMC)






