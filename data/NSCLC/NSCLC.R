source('/SAME/run_same.R')
source('/SAME/ComBat_sc.R')
source('/SAME/platform_batch_correction.R')

set.seed(2021)
NSCLC_seur <- readRDS("/data/NSCLC/NSCLC_seur.rds")
Idents(NSCLC_seur) <- "Celltype"
table(NSCLC_seur$Celltype)
NSCLC_seur <- RenameIdents(NSCLC_seur, `B cells` = "B", `Monocytes` = "Monocyte", `NK cells` = "NK cell", `NKT cells` = "T", `T cells CD4` = "T", `T cells CD8` = "T")
NSCLC_seur[["Celltype_used"]] <- Idents(NSCLC_seur)
#NSCLC_seur <- RenameAssays(NSCLC_seur, `RNA` = "TPM")
#NSCLC_seur[["RNA"]] <- CreateAssayObject((NSCLC_seur@assays$TPM@counts)/100) 
table(NSCLC_seur$Celltype_used)
saveRDS(NSCLC_seur, file = "data/NSCLC/NSCLC_seur.rds")

Y0 <- read.table("/data/NSCLC/Fig2b-WholeBlood_RNAseq.txt", row.names = 1, header = T, sep = "\t")
PBMC_seur <- readRDS("/data/NSCLC/PBMC_ImmuneCell_updated.rds")
Lung_seur <- readRDS("/data/NSCLC/Lung_ImmuneCell_updated.rds")

celltype.list <- c("B", "Monocyte", "NK cell", "T")
data_set_list <- list()
data_set_list[[1]] <- NSCLC_seur
data_set_list[[2]] <- subset(PBMC_seur, Celltype_used %in% celltype.list)
data_set_list[[3]] <- subset(Lung_seur, Celltype_used %in% celltype.list)
X <- data_set_list
Cl <- get_Cl(X)

YSG <- readRDS("/data/NSCLC/sg.list.rds")
tmp <- Reduce(intersect, lapply(X, function(x) rownames(x)))
YSG <- intersect(YSG, tmp)
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
W_tilde <- matrix(0.001, nrow = D, ncol = K)
seur_t0 <- data_set_list[[1]]
for (k in 1:K){
  seur_t0.list <- SplitObject(seur_t0, split.by = "Celltype_used")
  if (is.null(seur_t0.list[[celltype.list[k]]])){
    W_tilde[,k] = 0.001
  } else {
    seur_t0_k <- seur_t0.list[[celltype.list[k]]]
    W_tilde[,k] = apply(seur_t0_k@assays$RNA@counts[YSG,], 1, mean)
  }
}

mcmc_samples_theta1 = 100
Lambda = c(0:mcmc_samples_theta1)
rst <- Y_batch_correct(Y0, NSCLC_seur, YSG) # both Y0 and NSCLC_seur in TPM space
Y0_sc <- rst$bulk_to_sc
Y0_sc <- Y0_sc[YSG,] #TPM space
which(Y0_sc < 0)
tau_e_empirical <- rst$tau_e
Y_combind <- rst$before[YSG,]
#visualize the batch correction result
N = 12
pca_input <- cbind(rst$bulk_to_sc, rst$before[,13:24])
#pca_input <- as.matrix(rst$before)
batch <- c(rep(1, N), rep(2, N))
pca_input <- as.data.frame(t(pca_input))
pca_input <- scale(pca_input)
pca_rst<-prcomp(pca_input)
ggbiplot(pca_rst, obs.scale = 1, var.scale = 1,
         groups = as.factor(batch), ellipse = T, circle = T, var.axes = F) + 
  scale_color_discrete(name = '') +
  theme_bw()+ ggtitle("Y0 adjusted to single cell space")


W_tilde <- (W_tilde)/100
Y0_sc <- (Y0_sc)/100

#run same
rst_sc_alpha1 <- SAME(Y0_sc, X, W_tilde, mcmc_samples_theta1, Lambda, c_k, YSG, alpha =1)
rst_sc_alpha0 <- SAME(Y0_sc, X, W_tilde, mcmc_samples_theta1, Lambda, c_k, YSG, alpha =0)
rst_sc_alpha.5<- SAME(Y0_sc, X, W_tilde, mcmc_samples_theta1, Lambda, c_k, YSG, alpha =.5)
rst_y0_alpha.5<- SAME(Y0, X, W_tilde, mcmc_samples_theta1, Lambda, c_k, YSG, alpha =.5)

tmp <- t(rst_sc_alpha.5$theta1$z)
tmp_w <- rst_sc_alpha.5$theta1$w[[1]]
tmp2 <- t(rst_sc_alpha0$theta1$z)
tmp2_w <- rst_sc_alpha0$theta1$w[[1]]
tmp1 <- t(rst_y0_alpha.5$theta1$z)

colnames(tmp) <- celltype.list

WBC_gt <- read.table("/data/NSCLC/Fig2b_ground_truth_whole_blood.txt", row.names = 1, header = 1, sep = "\t")
gt_sub <- WBC_gt[,c(7,3,8,4)]
corr_tmp <- sapply(c(1:4), function(i) cor(gt_sub[,i], tmp[,i]))

cibersortx_s <- read.table("/data4/lbl//Yale/HCL/Bayesian/CibersortX/Fig2ab-NSCLC_PBMCs/CIBERSORTx_Adjusted_S.txt", header = T, row.names = 1, sep = "\t")
cibersortx_s$T <- rowSums(cibersortx_s[,c(1,3,4)])
cibersortx_s_sub <- cibersortx_s[,c(5,2,6,10)]
corr_cibersortx_s <- sapply(c(1:4), function(i) cor(gt_sub[,i], cibersortx_s_sub[,i]))
corr_df <- data.frame(corr = c(corr_tmp, corr_tmp2, corr_cibersortx_s), 
                      mode = c(rep("alpha0.5",4),rep("alpha0",4),rep("Smode",4)),
                      celltype = c(celltype.list, celltype.list, celltype.list))
p <- ggplot(corr_df, aes(x=mode, y=corr, colour = celltype)) + geom_jitter(width = 0.2)