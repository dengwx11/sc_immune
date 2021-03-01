library(rliger)
source("SAME/clean_format.R")

run_liger <- function(files,tissue_list,YSG,output_path){
    seur.TPM_list <- list()
    seur_list <- list()
    for(i in seq(length(tissue_list))){
        seur <- readRDS(files[i]) ## counts
        seur <- NormalizeData(seur,normalization.method = "LogNormalize",scale.factor = 1000000)
        seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^MT-")
        seur.TPM <- exp(GetAssayData(seur[['RNA']],slot='data'))-1
        seur.TPM_list[[i]] <- seur.TPM
        seur_list[[i]] <- seur
    }
    names(seur.TPM_list) <- tissue_list
    names(seur_list)<- tissue_list

    liger <- createLiger(seur.TPM_list) 
    liger <- normalize(liger)
    liger <- selectGenes(liger)
    liger@var.genes <- intersect(YSG,liger@var.genes)
    
    liger <- scaleNotCenter(liger)
    liger <- optimizeALS(liger, k = 20)
    liger <- quantile_norm(liger)
    liger <- louvainCluster(liger, resolution = 0.25)
    
    liger <- runTSNE(liger) 
    liger <- runUMAP(liger)
    
    saveRDS(liger,paste0(output_path,"/liger.rds"))
    
    ans<-list()
    ans$seur.TPM_list <- seur.TPM_list
    ans$seur_list <- seur_list
    ans$liger <- liger
    return(ans)
    
}

get_WH_W_tilde <- function(liger,seur.TPM_list,seur_list){
    WH <- list()
#    WVH <- list()
#    mat <- list()
    new_split <- list()
    tissue_list <- names(seur.TPM_list)

    for(i in seq(length(tissue_list))){
        WH[[i]] <- t(liger@W)%*%t(liger@H.norm[colnames(seur.TPM_list[[i]]),])
        #WH_PBMC <- t(ifnb_liger@W)%*%t(ifnb_liger@H$PBMC_HCL)
#        WVH[[i]]<- (t(liger@W)+t(liger@V[[1]]))%*%t(liger@H.norm[colnames(seur.TPM_list[[i]]),])
        #WVH_PBMC<- (t(ifnb_liger@W)+t(ifnb_liger@V$PBMC_HCL))%*%t(ifnb_liger@H$PBMC_HCL)
#        mat[[i]] <- seur.TPM_list[[i]][liger@var.genes,]
#        mat[[i]] <- as.matrix(mat[[i]])
        
        new_split[[i]] <- CreateSeuratObject(WH[[i]])
        new_split[[i]]$Celltype_used <- as.character(seur_list[[i]]$Celltype_used)
        new_split[[i]]$Batch <- rep(tissue_list[i],ncol(seur_list[[i]]))
    }
    
    W_tilde <- list()

    for(i in seq(length(tissue_list))){
        Idents(new_split[[i]]) <- "Celltype_used"
        new_split[[i]] <- NormalizeData(new_split[[i]])
        W_tilde[[i]] <- as.matrix(AverageExpression(new_split[[i]],slot='data')$RNA)
    }
    
    names(W_tilde) <- tissue_list
    
    return(W_tilde)

}


get_tissue_gene <- function(target_idx, ref_idx_list, W_tilde){
    ans <- list()
    for(i in seq(length(ref_idx_list))){
        ref_idx <- ref_idx_list[i]
        df.cor <- PlotCompare(target_idx,ref_idx,W_tilde)
        df.cor$ratio <- df.cor$tissue2/df.cor$tissue1

        df.cor.plot <- df.cor
        ratio_to_plot <- log(df.cor.plot$ratio)
        km <- Ckmeans.1d.dp(ratio_to_plot,k=3)

        df.cor.plot$cluster <- km$cluster
        ans[[i]] <- df.cor.plot[df.cor.plot$cluster %in% 2,]
    }
    return(ans)
}

