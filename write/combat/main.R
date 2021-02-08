#check NSCLC
sum(NSCLC_seur@assays$RNA@counts[,1])
cancer_mat <- as.matrix(NSCLC_seur@assays$RNA@counts)
sum(cancer_mat[,1])
N1 <- ncol(cancer_mat)

PBMC_ImmuneCell_updated <- readRDS("/data/PBMC_ImmuneCell_updated.rds")
#check PBMC and normalize raw counts to TPM space
PBMC_ImmuneCell_updated <- NormalizeData(PBMC_ImmuneCell_updated, scale.factor = 1000000)
sum(PBMC_ImmuneCell_updated@assays$RNA@counts[,1]) #4616 counts for the first cell
sum(expm1(PBMC_ImmuneCell_updated@assays$RNA@data[,1])) #1e+06
hcl_mat <- expm1(as.matrix(PBMC_ImmuneCell_updated@assays$RNA@data))
# PBMC_ImmuneCell_updated[["TPM"]] <- CreateAssayObject(counts = hcl_mat)
# PBMC_ImmuneCell_updated[["TPM"]] <- CreateAssayObject(data = log1p(hcl_mat))
sum(hcl_mat[,1])#1e+06
sum(PBMC_ImmuneCell_updated@assays$TPM@data[,1])
sum(PBMC_ImmuneCell_updated@assays$TPM@counts[,1])
N2 <- ncol(hcl_mat)

#run combat_sc
batch <- c(rep(1,N1), rep(2,N2))
gl <- intersect(rownames(cancer_mat), rownames(hcl_mat))

YSG <- intersect(YSG, gl)
mat <- cbind(cancer_mat[gl,], hcl_mat[gl,])
mat <- cbind(cancer_mat[YSG,], hcl_mat[YSG,])
rst <- ComBat_sc(log2(mat+1), batch)

################################
###### prepare SAME input ######
################################
celltype.list <- c("B", "Monocyte", "NK cell", "T")

Lung_ImmuneCell_updated <- NormalizeData(Lung_ImmuneCell_updated, scale.factor = 1000000)
hcl_mat <- expm1(as.matrix(Lung_ImmuneCell_updated@assays$RNA@data))
Lung_ImmuneCell_updated[["TPM"]] <- CreateAssayObject(counts = hcl_mat)
sum(Lung_ImmuneCell_updated@assays$TPM@counts[,1])

CB_ImmuneCell_updated <- NormalizeData(CB_ImmuneCell_updated, scale.factor = 1000000)
hcl_mat <- expm1(as.matrix(CB_ImmuneCell_updated@assays$RNA@data))
CB_ImmuneCell_updated[["TPM"]] <- CreateAssayObject(counts = hcl_mat)
sum(CB_ImmuneCell_updated@assays$TPM@counts[,1])

BM_ImmuneCell_updated <- NormalizeData(BM_ImmuneCell_updated, scale.factor = 1000000)
hcl_mat <- expm1(as.matrix(BM_ImmuneCell_updated@assays$RNA@data))
BM_ImmuneCell_updated[["TPM"]] <- CreateAssayObject(counts = hcl_mat)
sum(BM_ImmuneCell_updated@assays$TPM@counts[,1])

hcl_mat <- expm1(as.matrix(PBMC_ImmuneCell_updated@assays$RNA@data))
PBMC_ImmuneCell_updated[["TPM"]] <- CreateAssayObject(counts = hcl_mat)
sum(PBMC_ImmuneCell_updated@assays$TPM@counts[,1])
table(CB_ImmuneCell_updated$Celltype_used)
table(Liver_ImmuneCell_updated$Celltype_used) #no monocyte
table(Kidney_ImmuneCell_updated$Celltype_used) #no monocyte
table(BM_ImmuneCell_updated$Celltype_used)

NSCLC_seur[["TPM"]] <- CreateAssayObject(counts = 2^(rst$bulk_to_sc))
sum(NSCLC_seur@assays$TPM@counts[,1])
# #T0 = NSCLC, Tx = HCL
# data_set_list <- list()
# data_set_list[[1]] <- NSCLC_seur
# data_set_list[[2]] <- subset(PBMC_ImmuneCell_updated, Celltype_used %in% celltype.list)
# data_set_list[[3]] <- subset(CB_ImmuneCell_updated, Celltype_used %in% celltype.list)
# data_set_list[[4]] <- subset(BM_ImmuneCell_updated, Celltype_used %in% celltype.list)

#T0 = HCL_PBMC, Tx = HCL
# data_set_list <- list()
# data_set_list[[1]] <- subset(PBMC_ImmuneCell_updated, Celltype_used %in% celltype.list)
# data_set_list[[2]] <- subset(CB_ImmuneCell_updated, Celltype_used %in% celltype.list)
# data_set_list[[3]] <- subset(BM_ImmuneCell_updated, Celltype_used %in% celltype.list)


# #T0 = HCL_PBMC, Tx = NSCLC+HCL
data_set_list <- list()
data_set_list[[1]] <- subset(PBMC_ImmuneCell_updated, Celltype_used %in% celltype.list)
data_set_list[[2]] <- NSCLC_seur
data_set_list[[3]] <- subset(CB_ImmuneCell_updated, Celltype_used %in% celltype.list)
data_set_list[[4]] <- subset(BM_ImmuneCell_updated, Celltype_used %in% celltype.list)

X <- data_set_list
Cl <- get_Cl(X)
print(Cl$levels)
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
W_tilde <- matrix(0, nrow = D, ncol = K)
seur_t0 <- data_set_list[[1]]
for (k in 1:K){
  seur_t0.list <- SplitObject(seur_t0, split.by = "Celltype_used")
  if (is.null(seur_t0.list[[celltype.list[k]]])){
    W_tilde[,k] = 0
  } else {
    seur_t0_k <- seur_t0.list[[celltype.list[k]]]
    W_tilde[,k] = apply(seur_t0_k@assays$TPM@counts[YSG,], 1, mean)
  }
  colnames(W_tilde) <- celltype.list
}

mcmc_samples_theta1 = 50
Lambda = c(0:mcmc_samples_theta1)
#change line 13 code as mat <- X@assays$TPM@counts
Y0 <- read.table("/data/NSCLC/Fig2b-WholeBlood_RNAseq.txt", row.names = 1, header = T, sep = "\t")
# rst <- Y_batch_correct(Y0, NSCLC_seur, YSG) # both Y0 and NSCLC_seur in TPM space
rst <- Y_batch_correct(Y0, X[[1]], YSG)
Y0_sc <- rst$bulk_to_sc
Y0_sc <- Y0_sc[YSG,] #TPM space
N = ncol(Y0_sc)
C0 = sum(c_k)



rst_sc_alpha0 <- SAME(log1p(Y0_sc), X, log1p(W_tilde), mcmc_samples_theta1, Lambda, c_k, YSG, alpha =0)
rst_sc_alpha1 <- SAME(log1p(Y0_sc), X, log1p(W_tilde), mcmc_samples_theta1, Lambda, c_k, YSG, alpha =1)

# rst_sc_alpha.5<- SAME(log1p(Y0_sc), X, log1p(W_tilde), mcmc_samples_theta1, Lambda, c_k, YSG, alpha =.5)

#check readme.txt for details of the results
saveRDS(rst_sc_alpha1, file = "/wirte/combat/SAME_rst_t4_HCl_cancer_alpha1.rds")
saveRDS(rst_sc_alpha0, file = "/wirte/combat/SAME_rst_t3_HCl_cancer_alpha0.rds")

gamma_t2_1 <- Reduce("+", SAME_rst_t2_HCl_cancer_alpha1$theta2$gamma)
gamma_t2_2 <- Reduce("+", SAME_rst_t2_pi.8_HCl_cancer_alpha1$theta2$gamma)
gamma_t3_1 <- Reduce("+", SAME_rst_t3_allHCl_alpha1$theta2$gamma)
gamma_t3_2 <- Reduce("+", SAME_rst_t3_pi.8_allHCl_alpha1$theta2$gamma)

gamma_compare <- cbind(gamma_t2_1[,1], gamma_t2_2[,1], gamma_t3_1[,1], gamma_t3_2[,1])
gamma_compare <- data.frame(gamma_sum = c(gamma_t2_1[,1], gamma_t2_2[,1], gamma_t3_1[,1], gamma_t3_2[,1]),
                          para = c(rep("SAME_rst_t2_HCl_cancer_alpha1", 906), rep("SAME_rst_t2_pi.8_HCl_cancer_alpha1", 906), 
                                   rep("SAME_rst_t3_allHCl_alpha1", 906), rep("SAME_rst_t3_pi.8_allHCl_alpha1", 906)))

ggplot(data = gamma_compare, aes(x=gamma_sum))+
  geom_histogram(aes(fill=para))+
  facet_grid(.~para)