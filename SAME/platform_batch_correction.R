library(nnls)
library(sva)
set.seed(2020)
#generate pseudo bulk from single cell dataset, and return a result list of TPM pseudo bulk and cell type proprtion
generate_pseudobulk <- function(X, YSG){
  rst <- list()
  Cl <- levels(factor(as.vector(X$Celltype_used)))
  est_prop <- rep(0, length(Cl))
  meta_prop <- as.data.frame(table(X$Celltype_used))
  meta_prop[,2] <- meta_prop[,2]/sum(meta_prop[,2])
  rownames(meta_prop) <- meta_prop[,1]
  meta_prop <- meta_prop[,-1]
  mat <- X@assays$RNA_homolog@counts[YSG,]
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
  Y_pseudo <- (Y_pseudo/sum(Y_pseudo))*1000000
  rst$Y_pseudo <- Y_pseudo
  rst$prop <- est_prop
  return(rst)
}


#S-mode of CibersortX for single-cell batch correction; adjust the signature matrix
#This function returns adjusted W  
s_batch_correct <- function(Y0, w_t0, X){
  N <- ncol(Y0)
  Y_star <- sapply(c(1:N), function(i) generate_pseudobulk(X, YSG)$Y_pseudo)
  F_star <- t(sapply(c(1:N), function(i) generate_pseudobulk(X, YSG)$prop))
  batch <- c(rep(1, N), rep(2, N))
  Y_combine <- cbind(Y0, Y_star)
  Y_adj <- ComBat(dat = log2(Y_combine+1), batch = batch)
  Y_star_adj <- t(2^(Y_adj[,(N+1):(2*N)]))
  W_adj <- sapply(c(1:D), function(i)nnls(F_star, Y_star_adj[,i])$x)
  return(t(W_adj))
}

#B-mode of CibersortX for single-cell batch correction; adjust the bulk profile 
#This function returns adjusted Y0
b_batch_correct <- function(w_t0, Y0){
  N <- ncol(Y0)
  z_initial <- sapply(c(1:N), function(i) nnls(w_t0, Y0[,i])$x)
  Y <- w_t0 %*% z_initial
  Y_combine <- cbind(Y0, Y)
  batch <- c(rep(1, N), rep(2, N))
  Y_adj <- ComBat(dat = log2(Y_combine+1), batch = batch)
  Y0_adj <- (2^Y_adj[,1:N])-1
  return(Y0_adj)
}
