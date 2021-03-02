PBMC_seur <- readRDS("/data4/lbl/Yale/HCL/protein_coding/PBMC_ImmuneCell_updated.rds")
Lung_seur <- readRDS("/data4/lbl/Yale/HCL/protein_coding/Lung_ImmuneCell_updated.rds")
# Kidney_seur <- readRDS("/data4/lbl/Yale/HCL/protein_coding/Kidney_ImmuneCell_updated.rds")
BM_seur <- readRDS("/data4/lbl/Yale/HCL/protein_coding/BM_ImmuneCell_updated.rds")
CB_seur <- readRDS("/data4/lbl/Yale/HCL/protein_coding/CB_ImmuneCell_updated.rds")
# Liver_seur <- readRDS("/data4/lbl/Yale/HCL/protein_coding/Liver_ImmuneCell_updated.rds")

#signature gene, celltype, and tissue lists
YSG <- readRDS("/data4/lbl/Yale/HCL/SAME/sg.list.rds")
celltype.list <- c("B", "Neutrophil", "NK cell", "T")
tissue.list <- c("PBMC", "Lung", "CB", "BM")

#change expression profile into TPM space
# Lung_seur <- NormalizeData(Lung_seur, scale.factor = 1000000)
# hcl_mat <- expm1(as.matrix(Lung_seur@assays$RNA@data))
# Lung_seur[["TPM"]] <- CreateAssayObject(counts = hcl_mat)
# sum(Lung_seur@assays$TPM@counts[,1])

# CB_seur <- NormalizeData(CB_seur, scale.factor = 1000000)
# hcl_mat <- expm1(as.matrix(CB_seur@assays$RNA@data))
# CB_seur[["TPM"]] <- CreateAssayObject(counts = hcl_mat)
# sum(CB_seur@assays$TPM@counts[,1])

# BM_seur <- NormalizeData(BM_seur, scale.factor = 1000000)
# hcl_mat <- expm1(as.matrix(BM_seur@assays$RNA@data))
# BM_seur[["TPM"]] <- CreateAssayObject(counts = hcl_mat)
# sum(BM_seur@assays$TPM@counts[,1])

PBMC_seur <- NormalizeData(PBMC_seur, scale.factor = 1000000)
hcl_mat <- expm1(as.matrix(PBMC_seur@assays$RNA@data))
PBMC_seur[["TPM"]] <- CreateAssayObject(counts = hcl_mat)

############################
###prepare for SAME input###
############################
data_set_list <- list()
data_set_list[[1]] <- subset(PBMC_seur, Celltype_used %in% celltype.list)
data_set_list[[2]] <- subset(Lung_seur, Celltype_used %in% celltype.list)
data_set_list[[3]] <- subset(CB_seur, Celltype_used %in% celltype.list)
data_set_list[[4]] <- subset(BM_seur, Celltype_used %in% celltype.list)
rm(PBMC_seur, Lung_seur, CB_seur, BM_seur)

X <- data_set_list
Cl <- get_Cl(X)
print(Cl$levels)
tmp <- Reduce(intersect, lapply(X, function(x) rownames(x)))
Y0 <- read.table("/data/NSCLC/Fig2b-WholeBlood_RNAseq.txt", row.names = 1, header = T, sep = "\t")
YSG <- intersect(YSG, tmp)
YSG <- intersect(YSG, rownames(Y0))
T = length(data_set_list)
K = length(Cl$levels)
c_k = matrix(0, nrow = K, ncol = T)
D = length(YSG)

for (t in 1:T){
  seur = data_set_list[[t]]
  for (k in 1:K){
    seur.list <- SplitObject(seur, split.by = "Celltype_used")
    if (is.null(seur.list[[celltype.list[k]]])){
      c_k[k,t] = 0
    } else {
      c_k[k,t] = ncol(seur.list[[celltype.list[k]]])
    }
  }
}
W_tilde <- matrix(0, nrow = D, ncol = K)
seur_t0 <- data_set_list[[1]]
for (k in 1:K){
  seur_t0.list <- SplitObject(seur_t0, split.by = "Celltype_used")
  if (is.null(seur_t0.list[[celltype.list[k]]])){
    W_tilde[,k] = 0
  } else {
    seur_t0_k <- seur_t0.list[[celltype.list[k]]]
    W_tilde[,k] = apply(seur_t0_k@assays$RNA@counts[YSG,], 1, mean)
  }
  colnames(W_tilde) <- celltype.list
}

mcmc_samples_theta1 = 50
Lambda = c(0:mcmc_samples_theta1)
#change line 13 code as mat <- X@assays$TPM@counts

rst <- Y_batch_correct(Y0, X[[1]], YSG)
Y0_sc <- rst$bulk_to_sc
Y0_sc <- Y0_sc[YSG,] #TPM space
N = ncol(Y0_sc)
C0 = sum(c_k)

corr_celltype <- c()
pr_gamma_rst <- c()
empirical_pi <- c(0.1, 0.15, 0.2, 0.25, 0.3)
rst_new.list <- list()
for (i in 1:length(empirical_pi)) {
  rst_counts_alpha1 <- SAME(Y0_sc, X, W_tilde, empirical_pi[i], mcmc_samples_theta1, Lambda, c_k, YSG, alpha =1) #time  554.109s
  
  #calculate v*g
  est_vg <- lapply(c(1:mcmc_samples_theta1), function(i) rst_counts_alpha1$theta2$gamma[[i]] * rst_counts_alpha1$theta2$v[[i]])
  average_est_vg_1 <- Reduce('+',est_vg)/mcmc_samples_theta1
  rst_counts_alpha1$theta2$vg <- average_est_vg_1
  #estimate z by nnls
  z_est_vg_nnls <- t(sapply(c(1:N), function(j) nnls(average_est_vg_1,(Y0_sc[YSG,j]))$x))
  corr_vg_nnls<- sapply(c(1:4), function(i) cor(as.vector(as.matrix(gt_sub)[,i]), as.vector(z_est_vg_nnls[,i])))
  corr_vg_all <- cor(as.vector(as.matrix(gt_sub)), as.vector(z_est_vg_nnls))
  corr_all <- append(corr_all, corr_vg_all)
  #corr_celltype[i,] <- corr_vg_nnls
  rst_counts_alpha1$theta2$corr_celltype <- corr_vg_nnls
  rst_counts_alpha1$theta2$corr_all <- corr_vg_all
  rst_new.list[[i]] <- rst_counts_alpha1
  saveRDS(rst_counts_alpha1, file = paste0("SAME/V2/rst_tauv0.1_updatepi", empirical_pi[i], "_alpha1.rds"))
}
names(rst_new.list) <- paste0("pi",empirical_pi)

gamma_tmp <- lapply(1:5, function(j) unlist(lapply(rst_new.list[[j]]$theta2$gamma, function(x) sum(x))))
for (i in 1:5){
  corr_celltype <- c(corr_celltype,  rst_new.list[[i]]$theta2$corr_celltype)
  pr_gamma <- gamma_tmp[[i]]/(D*K)
  pr_gamma_rst <- c(pr_gamma_rst, pr_gamma)
}
corr_celltype <- c(corr_celltype, corr_wtilde)

df_corr_ct <- data.frame(corr = corr_celltype,
                         celltype = rep(celltype.list, 60),
                         initial_pi = rep(c(rep(rep(c(names(rst_new.list), "W_tilde"), each = 4), 5)), 2),
                         tau_v = rep(c(rep(0, 24), rep(1, 24), rep(0.1, 24), rep(4, 24), rep(0.01, 24)),2),
                         update_pi = c(rep("original", 120), rep("noupdate",120))
)
df_corr_ct$update_pi_tau_v <- paste0(df_corr_ct$update_pi, "_", df_corr_ct$tau_v)
saveRDS(df_corr_ct, file = "SAME/V2/corr_compare_df.rds")


p <- ggplot(df_corr_ct, aes(x=initial_pi, y=corr, colour = celltype)) + geom_jitter(width = 0.1) + facet_wrap(~update_pi_tau_v, nrow =2)

p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="red")
p1 <- ggplot(data=df_corr_ct,aes(x=initial_pi, y=corr, group=celltype, colour = celltype))+geom_line()

iteration = rep(c("1-10", "11-20", "21-30", "31-40", "41-50"), each = 10)
df_gamma <- data.frame(Pr = pr_gamma_rst,
                       initial_pi = rep(rep(rep(names(rst_new.list), each = 50), 5),2),
                       iteration = rep(rep(rep(iteration, 5),5),2),
                       tau_v = rep(c(rep(0, 250), rep(1, 250), rep(0.1, 250), rep(4, 250), rep(0.01, 250)),2),
                       update_pi = c(rep("original", 1250), rep("noupdate",1250))
)
df_gamma$update_pi_tau_v <- paste0(df_gamma$update_pi, "_", df_gamma$tau_v)
p2 <- ggplot(df_gamma, aes(x=initial_pi, y=Pr, colour = iteration)) + geom_jitter(width = 0.1)+ facet_wrap(~update_pi_tau_v, nrow = 2)
p2 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red")
saveRDS(df_gamma, file = "SAME/V2/gamma_compare_df.rds")
df_gamma <- readRDS(file = "SAME/V2/gamma_compare_df.rds")