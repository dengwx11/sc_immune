library(nnls)
library(MuSiC)
library(xbioc)
library(ggplot2)
library(Biobase)
library(plotROC)
library(Seurat)

run_music <- function(X_target,Y0,YSG, method="music", simulation=FALSE){
    ## X_target : Seurat Object
    
    D <- nrow(Y0)
    N <- ncol(Y0)
    
    assaycounts <- as.matrix(GetAssayData(X_target))[YSG,]
    pheno <- X_target@meta.data
    if(simulation){
        rownames(assaycounts) <- paste0('Gene',c(1:D))
        colnames(assaycounts) <- paste0('Cell',c(1:ncol(assaycounts)))
        pheno$SampleID <- paste0('Cell',c(1:ncol(assaycounts)))
        rownames(pheno) <- paste0('Cell',c(1:ncol(assaycounts)))
    }
    sc.eset <- ExpressionSet(assayData = assaycounts, phenoData = as(data.frame(pheno),"AnnotatedDataFrame"))
    
    
    assaycounts <- as.matrix(Y0)[YSG,]
    if(simulation){
        rownames(assaycounts) <- paste0('Gene',c(1:D))
        colnames(assaycounts) <- paste0('ID',c(1:N))
        pheno <- data.frame('subjectID' = paste0('ID',c(1:N)))
        rownames(pheno) <- paste0('ID',c(1:N))
    }else{
        pheno <- data.frame('subjectID' = colnames(Y0))
        rownames(pheno) <- colnames(Y0)
    }
    bulk.eset <- ExpressionSet(assayData = assaycounts, phenoData = as(data.frame(pheno),"AnnotatedDataFrame"))
    
    
    # Estimate cell type proportions
    Est.prop = music_prop(bulk.eset = bulk.eset, sc.eset = sc.eset, clusters = 'Celltype_used',
                                   samples = 'SampleID',  verbose = F)
    names(Est.prop)
    if(method=="music"){
        z_est <- t(Est.prop$Est.prop.weighted)
    }else if(method == "nnls"){
        z_est <- t(Est.prop$Est.prop.allgene)
    }
    
    
    return(z_est)
} 