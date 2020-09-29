get_Cl <- function(X){
    Cl <- Reduce(c, lapply(X, function(X) as.vector(X$Celltype_used)) )
    Cl <- factor(Cl)
    celltype_list <- levels(Cl)
    Cl <- as.numeric(Cl)
    Cl <- list(Cl = Cl, levels = celltype_list)
    return(Cl)
}

get_X_mat <- function(X, SG){
    X_mat <- list()
    for(i in 1:length(X)) X_mat[[i]] <- X[[i]]@assays$RNA@counts[YSG,]
    return(X_mat)
}