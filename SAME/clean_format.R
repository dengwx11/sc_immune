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


get_est_vg <- function(rst){
  est_vg <- lapply(c(1:mcmc_samples_theta1), function(i) rst$theta2$gamma[[i]] * rst$theta2$v[[i]])
  average_est_vg <- Reduce('+',est_vg)/mcmc_samples_theta1
  average_gamma_est <- Reduce('+',rst$theta2$gamma)/mcmc_samples_theta1
  ans <- list()
  ans$vg <-  average_est_vg
  ans$averg_gamma <- average_gamma_est
  return(ans)                 
}
                   
PlotCompare <- function(i,j,W_tilde){
    W_tilde_A <- W_tilde[[i]]
    W_tilde_B <- W_tilde[[j]]
    celltype.intersect <- intersect(colnames(W_tilde_A),colnames(W_tilde_B))
    df.cor1 <- melt(W_tilde_A[,celltype.intersect],value.name = "tissue1")
    df.cor2 <- melt(W_tilde_B[,celltype.intersect],value.name = "tissue2")
    df.cor <- cbind(df.cor1,tissue2=df.cor2[,3])
    colnames(df.cor)[c(1:2)] <- c("gene","celltype")
    return(df.cor)
}    
