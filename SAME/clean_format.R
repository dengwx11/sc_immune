get_used_cell <- function(celltype_lists,celltype_used_list){
    used_cell <- lapply(celltype_lists, function(celltype_list) which(celltype_list %in% celltype_used_list))
    return(used_cell) 
}
                        
                       

get_Cl <- function(celltype_lists,celltype_used_list){
    
    Cl <- Reduce(c, lapply(celltype_lists, function(celltype_list) as.vector(celltype_list)) )
    Cl <- factor(Cl,levels = celltype_used_list)
    Cl <- as.numeric(Cl)
    Cl <- list(Cl = Cl, levels = celltype_used_list)
    return(Cl)
}

get_X_mat <- function(X, YSG){
    
    X_mat <- lapply(X, function(x) x@assays$RNA@counts[YSG,])
    return(X_mat)
}

keep_used_cell <- function(seur_list,used_cell, YSG){
    seur_list <- lapply(seq(length(seur_list)), function(i) subset(seur_list[[i]], cells = Cells(seur_list[[i]])[used_cell[[i]]]))
    X_mat <- get_X_mat(seur_list, YSG)
    celltype_lists <- lapply(seq(length(seur_list)),function(i) seur_list[[i]]$Celltype_used)
    
    ans <- list()
    ans$X_mat <- X_mat
    ans$celltype_lists <- celltype_lists
    
    return(ans)                          
    
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
                   
get_c_k <- function(celltype_lists,celltype_used_list){
    
    c_k <- lapply(seq(length(celltype_lists)), function(i) table(celltype_lists[[i]]))
    c_k <- lapply(seq(length(celltype_lists)), function(i) c_k[[i]] <- c_k[[i]][celltype_used_list])             
    c_k <- t(Reduce(rbind,c_k))
    c_k[is.na(c_k)]<-0
    c_k <- t(c_k)
    rownames(c_k) <- names(celltype_lists)              
    return(c_k)              
}                  
                   
                   
