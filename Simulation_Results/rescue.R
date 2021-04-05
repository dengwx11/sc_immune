options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

T <- as.numeric(args[1])
D <- as.numeric(args[2])
K <- as.numeric(args[3])
corrupt_pi <- as.numeric(args[4])
tau_w_para <- as.numeric(args[5])
tau_xd_beta_para <- as.numeric(args[6])
mcmc_samples_theta1 <- as.numeric(args[7])
seed <- as.integer(args[8])
parameter <- as.character(args[9])


set.seed(seed)
setwd('/gpfs/ysm/pi/zhao-data/wd262/sc_immune/sc_immune')
source('Batch_Correction/platform_batch_correction_combat.R')
str_para = paste0("T=",T,".D=",D,".K=",K,".corrupt=",corrupt_pi,".tauW=",tau_w_para,".tauXdBeta=",tau_xd_beta_para)


if(parameter == "corrupt"){
    output_directory <- paste0("/gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/",parameter,"/",corrupt_pi*100)
}else if(parameter == "tissue_number"){
    output_directory <- paste0("/gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/",parameter,"/",T)
}else if(parameter == "gene_number"){
    output_directory <- paste0("/gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/",parameter,"/",D)
}else if(parameter == "unbalanced"){
    output_directory <- paste0("/gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/",parameter)
}
    
    
target_tissue <- '1'
tissue_list <- as.character(c(1:T))
celltype_used_list <- paste0("celltype ",c(1:K))
YSG <- as.character(c(1:D))
Y0 <- read.table(paste0(output_directory, "/Y0.txt"), sep = " ")
rownames(Y0) <- c(1:D)


mydir <- paste0("/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/simulation_rst/",parameter)            
rst <- readRDS(paste0(mydir,"/rst.",str_para,'.rds'))

rst_frac <- fraction_est(Y0, rst,target_tissue, celltype_used_list,method = c("tranSig","music","w_empirical"), adj.to.sc = F,simulation=T)


saveRDS(rst_frac,paste0(mydir,"/rst_frac.",str_para,'.rds'))