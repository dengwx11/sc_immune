# Update nuisance parameters, theta2
# output is matrix
set.seed(2020)



#update gamma
update_gamma <- function(W_T, pi_pre, v, tau_w, W_T_indicator
){
  #############
  #### Input:
  ## W_T: a list of W_hat from each tissue, and the length of the list is T, the dimension of each W_hat is D*K
  ## pi_pre: estimated pi from the last step. 
  ## v: D*K
  #### Output:
  ## gamma_new: a matrix of Bernoulli variables, D*K
  #############
  T = length(W_T)  
  gamma_new <- matrix(0, nrow = D, ncol = K)
  para <- log(pi_pre/(1-pi_pre)) - (tau_w/2)*Reduce("+", lapply(c(1:T), function(i)((W_T[[i]]-v)^2)*W_T_indicator[[i]])) + 
    (tau_w/2)*Reduce("+", lapply(c(1:T), function(x)W_T[[i]]^2*W_T_indicator[[i]]))
  for (d in 1:D){
    for (k in 1:K){
      p_temp <- 1/(1+exp(-para[d,k]))
      gamma_new[d,k] = rbinom(n=1, size = 1, prob = p_temp)
    }
  }
  return(gamma_new)
}

#update v
update_v <- function(tau_w, W_T, gamma,W_T_indicator,
                     tau_v=0.01
){
  #############
  #### Input:
  ## tau_w: variance of W_T, all w share the same tau_w
  ## W_tilde: empirical weight matrix, D*K
  ## W_T: a list of W_hat from each tissue, and the length of the list is T, the dimension of each W_hat is D*K
  ## gamma: a matrix of Bernoulli variables, D*K 
  #### Output:
  ## v_new: D*K
  #############                              
  v_new <- matrix(0, nrow = D, ncol = K)
  para1 <- T*tau_w + tau_v
  para2 <- (1/para1)*tau_w*Reduce("+", lapply(c(1:T), function(i) W_T[[i]]*W_T_indicator[[i]]))
  #print(c(para1,para2))
  #print(str(para2))
  for (d in 1:D){
    for (k in 1:K){
      if (gamma[d,k] == 1){
        v_new[d,k] <- rnorm(n=1, para2[d,k], sd = (1/sqrt(para1)))
      }
      if (gamma[d,k] == 0){
        v_new[d,k] <- rnorm(n=1, 0, sd = 1/sqrt(tau_v))
      }
    }
  }
  return(v_new)
}

#update pi_est
update_pi_est <- function(gamma, 
                          alpha_pi = 10^(-6), beta_pi = 10^(-6)
){                                                                
  para1 <- alpha_pi + sum(gamma)
  para2 <- beta_pi + D*K - sum(gamma)
  pi_new <- rbeta(n=1, para1, para2)
  return(pi_new)
}

#update tau_x
update_tau_x <- function(X, #X is a list of gene expression matrix from each tissue. For each matrix, cells should be ordered by cell types
                         W_T_indicator,
                         W_T, Cl, c_k,
                         alpha_x = 10^(-6), beta_x = 10^(-6)
){
  #############
  #### Input:
  ## X: a list of scRNA-seq gene expression matrix from each tissue. the dimension of X[[t]] is D*sum(c_k[,t])
  ## W_T: a list of W_hat from each tissue, and the length of the list is T, the dimension of each W_hat is D*K
  ## c_k: sequence of cell counts for each cell type and each tissue, K*T
  #### Output:
  ## tau_x_new: 1*D
  #############
  C0 <- Reduce('+',lapply(W_T_indicator,function(x) sum(x)))
    
  tau_x_new <- matrix(0,nrow = 1, ncol = D)
  para2_list <- list()
  para2_temp <- matrix(0, nrow = 1, ncol = D)
  end_idx = 0
  for (t in 1:T){
    w_t <- W_T[[t]]
    w_t_indicator <- W_T_indicator[[T]]
    w_temp <- matrix(0, nrow = D, ncol = sum(c_k[,t]))
    w_indicator_temp <- matrix(0, nrow = D, ncol = sum(c_k[,t]))
    start_idx <- end_idx + 1
    end_idx <- start_idx + sum(c_k[,t])-1
    for (k in 1:K){                         
      idx = which(Cl[start_idx:end_idx] == k)
      w_temp[,idx] = w_t[,k]
      w_indicator_temp[,idx] = w_t_indicator[,k]
    }
    x_t <- as.matrix(X[[t]])
    para2_term1 <- ((x_t - w_temp)^2*w_indicator_temp)/2
    para2_list[[t]] <- apply(para2_term1, 1, sum)
    para2_temp = para2_temp + para2_list[[t]]
  }
  para1 <- rep(((C0/2) + alpha_x), D)
  para2 <- para2_temp + beta_x
  for(d in 1:D){
    tau_x_new[,d] <- rgamma(1, para1[d],para2[d])
  }
  return(tau_x_new)
}

#update tau_w
update_tau_w <- function (W_T, v, gamma, W_T_indicator,
                          alpha_w = 10^(-6), beta_w = 10^(-6)
){
  #############
  #### Input:
  ## W_T: a list of W_hat from each tissue, and the length of the list is T, the dimension of each W_hat is D*K
  ## v: D*K
  ## gamma: a matrix of Bernoulli variables, D*K 
  #### Output:
  ## tau_w_new: all the w share this same tau_w
  #############                                   
  para1 <- (T*D*K/2) + alpha_w
  para2 <- sum(Reduce("+", lapply(c(1:T), function(i)(W_T[[i]]-v*gamma)^2*W_T_indicator[[i]])))/2 + beta_w
  tau_w_new <- rgamma(1, para1, para2)
  return(tau_w_new)
}



# library(Seurat)
# #preprocess scRNA-seq data to get X
# data_preprocessing <- function(data_set_list, # a list of scRNA-seq dataset, and the first item is t0 tissue. All the items are Seurat object.
#                                 CL = "Celltype_used", #for HCL data
#                                 SG, # a list of signature genes
#                                 celltype.list = c("B", "Dendritic cell", "Macrophage", "Monocyte", "Neutrophil", "NK cell", "Plasmocyte", "T")

# ){
#     X_SAME <- list()
#     T = length(data_set_list)
#     K = length(celltype.list)
#     c_k = matrix(0, nrow = K, ncol = T)
#     D = length(SG)
#     for (t in 1:T){
#         seur = data_set_list[[t]]
#         mat = matrix(0, nrow = D)
#         for (k in 1:K){
#             seur.list <- SplitObject(seur, split.by = CL)
#             if (is.null(seur.list[[celltype.list[k]]])){
#                 c_k[k,t] = 0
#             } else {
#                 c_k[k,t] = ncol(seur.list[[celltype.list[k]]])
#                 seur_k <- seur.list[[celltype.list[k]]]
#                 mat <- cbind(mat, seur_k@assays$RNA@counts[SG,])
#             }
#         }
#         X_SAME[[t]] = mat
#     }
#     W_tilde <- matrix(0, nrow = D, ncol = K)
#     seur_t0 <- data_set_list[[1]]
#     for (k in 1:K){
#         seur_t0.list <- SplitObject(seur_t0, split.by = CL)
#         if (is.null(seur.list[[celltype.list[k]]])){
#             W_tilde[,k] = 0
#         } else {
#             seur_t0_k <- seur_t0.list[[celltype.list[k]]]
#             W_tilde[,k] = apply(seur_k@assays$RNA@counts[SG,], 1, mean)
#         }
#     }
#     rst <- list()
#     rst$X <- X_SAME
#     rst$T <- T
#     rst$K <- K
#     rst$D <- D
#     rst$c_k <- c_k
#     rst$W_tilde <- W_tilde
#     return(rst)
# }

# #empirical gamma based on W_tilde
# update_gamma_empirical <- function(W_tilde, 
#                          W_T, pi_pre, v, tau_w
# ){
#     #############
#     #### Input:
#     ## W_tilde: empirical weight matrix, D*K
#     ## W_T: a list of W_hat from each tissue, and the length of the list is T, the dimension of each W_hat is D*K
#     ## pi_pre: estimated pi from the last step. 
#     ## v: D*K
#     #### Output:
#     ## gamma_new: a matrix of Bernoulli variables, D*K
#     #############                                 
#     gamma_new <- matrix(1, nrow = D, ncol = K)
#     para <- log(pi_pre/(1-pi_pre)) - (tau_w/2)*Reduce("+", lapply(W_T, function(x)(x-v)^2)) + 
#         (tau_w/2)*Reduce("+", lapply(W_T, function(x)x^2))
#     for (d in 1:D){
#         for (k in 1:K){
#             if(W_tilde[d,k] == 0){
#                 p_temp <- 1/(1+exp(-para[d,k]))
#                 gamma_new[d,k] = rbinom(n=1, size = 1, prob = p_temp)
#             }
#         }
#     }
#     return(gamma_new)
# }

