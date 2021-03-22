source("SAME/ComBat_sc.R")
library(nnls)
library(sva)
set.seed(2020)
#generate pseudo bulk from single cell dataset, and return a result list of TPM pseudo bulk and cell type proprtion
generate_pseudobulk <- function(X, YSG, Celltype_used = "Celltype_used"){
  rst <- list()
  Cl <- levels(factor(as.vector(X$Celltype_used)))
  est_prop <- rep(0, length(Cl))
  meta_prop <- as.data.frame(table(X$Celltype_used))
  meta_prop[,2] <- meta_prop[,2]/sum(meta_prop[,2])
  rownames(meta_prop) <- meta_prop[,1]
  meta_prop <- meta_prop[,-1]
  mat <- exp(X@assays$RNA@data[YSG,])-1
  #psb_mat <- matrix(0, ncol = 10001, nrow = length(YSG))
  psb_mat <- data.frame(rep(0, length(YSG)), row.names = YSG)
  for(i in 1:length(Cl)){
    tmp <- rnorm(1, mean = meta_prop[i], sd = 2*meta_prop[i])
    if(tmp > 0){
      est_prop[i] = tmp
      cells <- sample(rownames(as.data.frame(X$Celltype_used[X$Celltype_used == Cl[i]])), round(10000*tmp), replace = T)
      psb_mat <- cbind(psb_mat, as.data.frame(mat[,cells]))}
  }
  est_prop <- est_prop/sum(est_prop)
  Y_pseudo <- apply(psb_mat, 1, sum)
  #Y_pseudo <- (Y_pseudo/sum(Y_pseudo))*1000000
  Y_pseudo <- Y_pseudo/10000 #if X in TPM space
  rst$Y_pseudo <- Y_pseudo
  rst$prop <- est_prop
  return(rst)
}


#S-mode of CibersortX for single-cell batch correction; adjust the signature matrix
#This function returns adjusted W  
s_batch_correct <- function(Y0, w_t0, X, YSG, Celltype_used = "Celltype_used"){
  N <- ncol(Y0)
  rst.list <- sapply(c(1:N), function(i) generate_pseudobulk(X, YSG))
  Cl <- levels(factor(as.vector(X[[Celltype_used]])))
  Y_idx <- seq(1, 2*N, 2)
  F_idx <- Y_idx + 1
  Y.list <- rst.list[Y_idx]
  F.list <- rst.list[F_idx]
  Y_star <- Reduce(cbind, Y.list)
  F_star <- t(Reduce(cbind, F.list))
  batch <- c(rep(1, N), rep(2, N))
  gene.list <- intersect(rownames(Y0), YSG)
  Y_combine <- cbind(Y0[gene.list,], Y_star[gene.list,])
  Y_adj <- ComBat(dat = log2(Y_combine+1), batch = batch
                    #par.prior = FALSE, mean.only = TRUE
                    )
  Y_star_adj <- t(2^(Y_adj[,(N+1):(2*N)])-1)
  W_adj <- sapply(c(1:D), function(i)nnls(F_star, Y_star_adj[,i])$x)
  W_adj = t(W_adj)
  rownames(W_adj) <- gene.list
  colnames(W_adj) <- Cl
  return(t(W_adj))
}

#adjust Y0 by pseudo bulk from single cell reference
Y0_batch_correct <- function(Y0, X, YSG, Celltype_used = "Celltype_used"){
  N <- ncol(Y0)
  rst.list <- sapply(c(1:N), function(i) generate_pseudobulk(X, YSG))
  Cl <- levels(factor(as.vector(X$Celltype_used)))
  Y_idx <- seq(1, 2*N, 2)
  F_idx <- Y_idx + 1
  Y.list <- rst.list[Y_idx]
  F.list <- rst.list[F_idx]
  Y_star <- Reduce(cbind, Y.list)
  F_star <- t(Reduce(cbind, F.list))
  batch <- c(rep(1, N), rep(2, N))
  gene.list <- intersect(rownames(Y0), YSG)
  Y_combine <- cbind(Y0[gene.list,], Y_star[gene.list,])
  Y_adj <- ComBat(dat = log2(Y_combine+1), batch = batch
                    #par.prior = FALSE, mean.only = TRUE
                    )
#   Y_star_adj <- t(2^(Y_adj[,(N+1):(2*N)])-1)
#   W_adj <- sapply(c(1:D), function(i)nnls(F_star, Y_star_adj[,i])$x)
#   W_adj = t(W_adj)
#   rownames(W_adj) <- gene.list
#   colnames(W_adj) <- Cl
  Y0_adj <- 2^(Y_adj[,1:N])-1
  # rst <- list()
  # rst$before <- Y_combine 
  # rst$after <- Y_adj #Y0_adj = Y_adj[1:N, ]
  # return(rst)
  return(Y0_adj)
}

#B-mode of CibersortX for single-cell batch correction; adjust the bulk profile 
#This function returns adjusted Y0
b_batch_correct <- function(w_t0, Y0){
  N <- ncol(Y0)
  gene.list <- intersect(rownames(w_t0), rownames(Y0))
  z_initial <- sapply(c(1:N), function(i) nnls(w_t0[gene.list,], Y0[gene.list,i])$x)
  Y <- w_t0[gene.list,] %*% z_initial
  Y_combine <- cbind(Y0[gene.list,], Y)
  batch <- c(rep(1, N), rep(2, N))
  Y_adj <- ComBat(dat = log2(Y_combine+1), batch = batch)
  Y0_adj <- (2^Y_adj[,1:N])-1
  return(Y0_adj)
}


#adjust Y0 by pseudo bulk from single cell reference
Y_batch_correct <- function(Y0, X, YSG, Celltype_used = "Celltype_used"){
  N <- ncol(Y0)
  rst.list <- sapply(c(1:N), function(i) generate_pseudobulk(X, YSG))
  Cl <- levels(factor(as.vector(X[[Celltype_used]])))
  Y_idx <- seq(1, 2*N, 2)
  F_idx <- Y_idx + 1
  Y.list <- rst.list[Y_idx]
  F.list <- rst.list[F_idx]
  Y_star <- Reduce(cbind, Y.list)
  F_star <- t(Reduce(cbind, F.list))
  batch <- c(rep(1, N), rep(2, N))
  gene.list <- intersect(rownames(Y0), YSG)
  Y_combine <- cbind(Y0[gene.list,], Y_star[gene.list,])
  # Y_adj <- ComBat(dat = log2(Y_combine+1), batch = batch)
  # rst <- list()
  # rst$before <- Y_combine 
  # rst$after <- (2^Y_adj)-1 #Y0_adj = Y_adj[1:N, ]
  rst <- ComBat_sc(dat = as.matrix(log2(Y_combine+1)), batch = batch)
  rst$bulk_to_sc <- (2^rst$bulk_to_sc)
  rst$before <- Y_combine
  rst$YSG <- gene.list
  return(rst)
}

#estimate cell type fractions by nnls.
fraction_est <- function(Y0, SAME_rst, 
                        target_tissue, celltype_used_list,Celltype_used = "Celltype_used",
                        method = "nnls", adj.to.sc = TRUE
){
  YSG <- SAME_rst$YSG
  vg <- SAME_rst$vg
  rownames(vg) <- YSG
  X <- SAME_rst$X[[target_tissue]]
  N <- ncol(Y0)
  if(adj.to.sc){
    rst_bulk <- Y_batch_correct(Y0, X, YSG)
    YSG <- rst_bulk$YSG
    Y0_adj <- rst_bulk$bulk_to_sc
    Y0_adj <- Y0_adj[YSG,] #TPM space
  }else{
    YSG <- intersect(YSG, rownames(Y0))
    Y0_adj <- Y0[YSG,]
  }
  if(method == "nnls"){
    z_est <- t(sapply(c(1:N), function(n) nnls(vg[YSG,],(Y0_sc[YSG,n]))$x))
    colnames(z_est) <- celltype_used_list
  }
  rst <- list()
  rst$z_est <- z_est
  rst$YSG <- YSG
  return(rst)
}