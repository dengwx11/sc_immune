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
#parameter <- "corrupt" or "tissue_number" or "gene_number"
#T=3
#D=500
#K=8
#corrupt_pi <-.25
#tau_w_para <-1
#tau_xd_beta_para <-1
#mcmc_samples_theta1 <-40
#seed <-101
#parameter <-'tissue_number'

print(mcmc_samples_theta1)

set.seed(seed)
setwd('/gpfs/ysm/pi/zhao-data/wd262/sc_immune/sc_immune')
source('SAME/run_same.v2.R')
source('SAME/run_Simulation.R')
source('Simulation_Results/prepare_simulation.R')
source('Batch_Correction/platform_batch_correction_combat.R')

pi_ber = 0.3
N = 200 # bulk Y sample size
Iteration = 500 ## iteration number to get the largest angle between the vectors
str_para = paste0("T=",T,".D=",D,".K=",K,".corrupt=",corrupt_pi,".tauW=",tau_w_para,".tauXdBeta=",tau_xd_beta_para,".seed=",seed)

same_input <- generate_same_input(T,D,K,pi_ber,N,Iteration,corrupt_pi=corrupt_pi, unbalanced = TRUE)
if(parameter == "corrupt"){
    output_directory <- paste0("/gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/",parameter,"/",corrupt_pi*100)
}else if(parameter == "tissue_number"){
    output_directory <- paste0("/gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/",parameter,"/",T)
}else if(parameter == "gene_number"){
    output_directory <- paste0("/gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/",parameter,"/",D)
}else if(parameter == "unbalanced"){
    output_directory <- paste0("/gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/",parameter,"/",corrupt_pi*100)
}
    
dir.create(output_directory)    
output_directory <- paste0(output_directory,"/",seed)
dir.create(output_directory)
prepare_same_input_in_simulation(same_input, output_directory)

target_tissue <- '1'
files = list.files(path = output_directory,pattern = "X_*", full.names = TRUE)
tissue_list <- as.character(c(1:T))
celltype_used_list <- paste0("celltype ",c(1:K))
YSG <- as.character(c(1:D))
output_path = ""
Y0 <- read.table(paste0(output_directory, "/Y0.txt"), sep = " ")
rownames(Y0) <- c(1:D)


rst <- run_SAME(target_tissue,tissue_list,celltype_used_list,files,YSG,empirical_pi=0.3,mcmc_samples_theta1,
                liger.turnon=FALSE,output_path)

rst_frac <- fraction_est(Y0, rst,target_tissue, celltype_used_list,method = c("tranSig","music","w_empirical"), adj.to.sc = F,simulation=T)

rst[['X']] <- NULL
mydir <- paste0("/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/simulation_rst/",parameter)            
saveRDS(rst,paste0(mydir,"/rst.",str_para,'.rds'))
saveRDS(rst_frac,paste0(mydir,"/rst_frac.",str_para,'.rds'))