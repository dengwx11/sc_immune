library(Seurat)
source('SAME/run_same.R')
source("SAME/Simulation.R")
#preprocess scRNA-seq data to get X
data_preprocessing <- function(data_set_list, # a list of scRNA-seq dataset, and the first item is t0 tissue. All the items are Seurat object.
                                CL = "SAME_celltype", #for HCL data
                                SG, # a list of signature genes
                                celltype.list = c("B", "Dendritic cell", "Macrophage", "Monocyte", "Neutrophil", "NK cell", "Plasmocyte", "T")

){
    X_SAME <- list()
    T = length(data_set_list)
    K = length(celltype.list)
    c_k = matrix(0, nrow = K, ncol = T)
    D = length(SG)
    for (t in 1:T){
        seur = data_set_list[[t]]
        mat = matrix(0, nrow = D)
        for (k in 1:K){
            seur.list <- SplitObject(seur, split.by = CL)
            if (is.null(seur.list[[celltype.list[k]]])){
                c_k[k,t] = 0
            } else {
                c_k[k,t] = ncol(seur.list[[celltype.list[k]]])
                seur_k <- seur.list[[celltype.list[k]]]
                mat <- cbind(mat, seur_k@assays$RNA@counts[SG,])
            }
        }
        X_SAME[[t]] = mat
    }
    W_tilde <- matrix(0.001, nrow = D, ncol = K)
    seur_t0 <- data_set_list[[1]]
    for (k in 1:K){
        seur_t0.list <- SplitObject(seur_t0, split.by = CL)
        if (is.null(seur_t0.list[[celltype.list[k]]])){
            W_tilde[,k] = 0.001
        } else {
            seur_t0_k <- seur_t0.list[[celltype.list[k]]]
            W_tilde[,k] = apply(seur_t0_k@assays$RNA@counts[SG,], 1, mean)
        }
    }
    rst <- list()
    rst$X <- X_SAME
    rst$T <- T
    rst$K <- K
    rst$D <- D
    rst$c_k <- c_k
    rst$W_tilde <- W_tilde
    return(rst)
}
##load dataset 
data_dir <- "~/Mydata/HCL/SAME/"
data_set_list <- list()
filenames <- c("Lung", "PBMC", "BM", "Liver", "CB", "Kidney")
for (i in 1:length(filenames)){
  file_dir <- paste0(data_dir, filenames[i], "_sg_seur.rds")
  seur <- readRDS(file_dir)
  data_set_list[[i]] <- seur
}
##subset seur object with 3 cell types
data_set_list_CT3 <- list()
filenames <- c("Lung", "PBMC", "BM", "Liver", "CB", "Kidney")
for (i in 1:length(filenames)){
  file_dir <- paste0(data_dir, filenames[i], "_sg_seur.rds")
  seur <- readRDS(file_dir)
  seur <- subset(seur, subset = SAME_celltype == "Macrophage"| SAME_celltype == "Dendritic"|SAME_celltype == "T")
  data_set_list_CT3[[i]] <- seur
}

#data_set_list_CT5 <- list()
#filenames <- c("Lung", "PBMC", "BM", "Liver", "CB", "Kidney")
#for (i in 1:length(filenames)){
#  file_dir <- paste0(data_dir, filenames[i], "_sg_seur.rds")
#  seur <- readRDS(file_dir)
#  seur <- subset(seur, subset = SAME_celltype == "Macrophage"| SAME_celltype == "Dendritic"|SAME_celltype == "T" |SAME_celltype == "Neutrophil"|SAME_celltype == "Plasmocyte")
#  data_set_list_CT5[[i]] <- seur
#}
#PseudoBulk_CT3 <- list()
#Proportion_CT3 <- list()
#PseudoBulk_CT5 <- list()

##Lung pseudo bulk
PseudoBulk <- data.frame()
PseudoBulk_var <- data.frame()
seur <- data_set_list_CT3[[1]]
mat <- seur@assays$RNA@counts
meta_df <- seur[["SAME_celltype"]]
Lung_stat <- as.data.frame(table(meta_df))
Lung_stat[,2] = Lung_stat[,2]/sum(Lung_stat[,2])
set.seed(2020)
times = 50
for(i in 1:times){
  cells <- sample(Cells(seur), 10000, replace = T)
  mat_new <- mat[,cells]
  meta_df_new <- meta_df[cells,]
  stat <- as.data.frame(table(meta_df_new))
  Lung_stat[,i+2] =stat[,2]/10000
  for (j in 1:length(rownames(mat_new))) {
    PseudoBulk[j,i] = sum(mat_new[j,]) 
    PseudoBulk_var[j,i] =  var(mat_new[j,])
  }
  print(paste(i, "th times"))
}

rownames(Lung_stat) = Lung_stat[,1]
Lung_stat = Lung_stat[,-1]

#for (i in 1:6){
#  PseudoBulk <- data.frame()
#  seur <- data_set_list_CT3[[i]]
#  mat <- seur@assays$RNA@counts
#  meta_df <- seur[["SAME_celltype"]]
#  Lung_stat <- as.data.frame(table(meta_df))
#  Lung_stat[,2] = Lung_stat[,2]/sum(Lung_stat[,2])
#  set.seed(2020)
#  times = 100
#  for(i in 1:times){
#    cells <- sample(Cells(seur), 10000, replace = T)
#    mat_new <- mat[,cells]
#    meta_df_new <- meta_df[cells,]
#    stat <- as.data.frame(table(meta_df_new))
#    Lung_stat[,i+2] =stat[,2]/10000
#    for (j in 1:length(rownames(mat_new))) {
#      PseudoBulk[j,i] = sum(mat_new[j,]) + var(mat_new[j,])/50
#    }
#    #print(paste(i, "th times"))
#  }
#  rownames(Lung_stat) = Lung_stat[,1]
#  Lung_stat = Lung_stat[,-1]
#  PseudoBulk_CT3[[i]] <- PseudoBulk
#  Proportion_CT3[[i]] <- Lung_stat
#}
####SAME
sg.list <- readRDS("~/Mydata/HCL/SAME/sg.list.rds")
X <- data_set_list_CT3 ##change the dataset here
Cl <- get_Cl(X)
SAME_Input <- data_preprocessing(data_set_list = data_set_list_CT3, ##change the dataset here
                                  CL = "SAME_celltype",SG = sg.list, celltype.list = Cl$levels)

W_tilde <- SAME_Input$W_tilde
mcmc_samples_theta1 = 100
Lambda = c(0:mcmc_samples_theta1)
c_k <- SAME_Input$c_k
YSG <- sg.list
T <- SAME_Input$T
K <- SAME_Input$K
D <- SAME_Input$D
#Y0 <- readRDS("~/Mydata/HCL/SAME/PseudoBulk_Lung_sg_CT3.rds")
#Y0 <- as.matrix(Y0)/10000 #divided by the number of cells
Pseudo_count <- readRDS("~/Mydata/HCL/SAME/PseudoBulk_Lung_sg_CT3_S50.rds")
Pseudo_var <- readRDS("~/Mydata/HCL/SAME/PseudoBulk_Lung_sg_CT3_S50_var.rds")
N = ncol(Pseudo_count)
Pseudo_error_50 <- matrix(0, nrow = D, ncol = N)
for (d in 1:D){
  for (n in 1:N){
    Pseudo_error_50[d,n] <- rnorm(1, mean = 0, sd = sqrt(Pseudo_var[d,n]/50))
  }
}
Pseudo_norm_var50 <- (as.matrix(Pseudo_count)/10000) + Pseudo_error_50
Pseudo_norm_var50[Pseudo_norm_var50 <= 0.001 ] = 0.001
Y0_var50 <- Pseudo_norm_var50


##run SAME 
rst_var50_alpha1 <- SAME(Y0_var50, X, W_tilde, mcmc_samples_theta1, Lambda, c_k, YSG, alpha =1)
rst_var50_alpha0 <- SAME(Y0_var50, X, W_tilde, mcmc_samples_theta1, Lambda, c_k, YSG, alpha =0)
rst_var50_alpha0.5 <- SAME(Y0_var50, X, W_tilde, mcmc_samples_theta1, Lambda, c_k, YSG, alpha =0.5)
par(mfrow = c(1,2))
plot(as.vector(rst_var50_alpha1$theta1$z), as.matrix(Lung_sg_CT5_S50_stat[,-1]), xlab = "z_lambda100", ylab = "true_z", main = "alpha = 1, e_var = var/50")
plot(as.vector(rst_var50_alpha1$theta1$w[[1]]), W_tilde, xlab = "w_t0_lambda100", ylab = "w_tilde", main = "alpha = 1, e_var = var/50")


plot(as.vector(rst_var50_alpha0$theta1$z), as.matrix(Lung_sg_CT5_S50_stat[,-1]), xlab = "z_lambda100", ylab = "true_z", main = "alpha = 0, e_var = var/50")
plot(as.vector(rst_var50_alpha0$theta1$w[[1]]), W_tilde, xlab = "w_t0_lambda100", ylab = "w_tilde", main = "alpha = 0, e_var = var/50")


plot(as.vector(rst_var50_alpha0.5$theta1$z), as.matrix(Lung_sg_CT5_S50_stat[,-1]), xlab = "z_lambda100", ylab = "true_z", main = "alpha = 0.5, e_var = var/50")
plot(as.vector(rst_var50_alpha0.5$theta1$w[[1]]), W_tilde, xlab = "w_t0_lambda100", ylab = "w_tilde", main = "alpha = 0.5, e_var = var/50")





##generate pseudo bulk 
alpha = 0.5
w_pseudo <- alpha * W_tilde + (1-alpha) * w_t0_100
tau_e = 1000
mu_e = 0
Z_sim <- generate_z(K, N)
e <- matrix(rnorm(D*N, mu_e, sd = (1/sqrt(tau_e))), nrow = D, ncol = N)
e[e < 0] = 0.001
Y_sim = w_pseudo %*% Z_sim + e
saveRDS(Y_sim, "~/Mydata/HCL/SAME/PseudoBulk_Lung_sg_CT3_0.5w_t0.rds")
saveRDS(Z_sim, "~/Mydata/HCL/SAME/Lung_sg_CT5_0.5w_t0_stat.rds")
saveRDS(w_t0_100, "~/Mydata/HCL/SAME/Lung_w_t0_lambda100_CT3.rds")

coherent.genes <- c("CD68", "CD163", "CD3E", "CD4", "CD8A", "FOXP3", "GZMA", "NKG7", "CD19", "CD79A", "CD79B", "CD38", "FCER2", "MZB1", "FCER1A", "CPA3", "CD33", "CEACAM8", "CD209", "LYZ", "CD14", "FCGR3A")
pdf("coherent.genes.pdf", width = 15)

for (i in 1:length(coherent.genes)){
  for (j in 1:6){
    assign(paste0("p", j), VlnPlot(data_set_list[[j]], features = coherent.genes[i], group.by = "SAME_celltype", pt.size = 0))
  }
  print(plot_grid(p1,p2,p3,p4,p5,p6, ncol = 3))
}
dev.off()





