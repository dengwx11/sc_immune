library(rliger)
library(Seurat)
set.seed(2021)

## Data Loading
NSCLC_seur <- readRDS("data/NSCLC/NSCLC_seur.rds")



PBMC_seur <- readRDS("data/NSCLC/PBMC_ImmuneCell_updated.rds")
Lung_seur <- readRDS("data/NSCLC/Lung_ImmuneCell_updated.rds")

## Stage I: preprocessing and normalization
ifnb_liger <- createLiger(list(NSCLC = GetAssayData(NSCLC_seur[['RNA']],slot='counts'), 
                               PBMC_HCL = GetAssayData(PBMC_seur[['RNA']],slot='counts')))
ifnb_liger <- normalize(ifnb_liger)
ifnb_liger <- selectGenes(ifnb_liger) 
ifnb_liger <- scaleNotCenter(ifnb_liger) 

## Stage II: joint matrix factorization
ifnb_liger <- optimizeALS(ifnb_liger, k = 20)

## Stage III: quantile normalization and joint clustering