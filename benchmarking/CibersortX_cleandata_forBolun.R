## tissue_number
files_input = list.files(path = '/gpfs/ysm/scratch60/cpsc424/wd262/sc_immune/simulation/tissue_number',pattern = "*", full.names = TRUE)
tissue_number_list <- sapply(files_input, basename)
files_input <- as.vector(sapply(files_input, function(file) paste0(file,'/',c(101:115),'/')))
files_W_empirical = as.vector(sapply(files_input, function(file) paste0(file,'W_empirical.txt')))
files_W_transig = as.vector(sapply(files_input, function(file) paste0(file,'W_transig.txt')))
files_Y0 = as.vector(sapply(files_input, function(file) paste0(file,'Y0.txt')))                                    


names(tissue_number_list) <- NULL
tissue_number_list <- as.numeric(tissue_number_list)
#tissue_number_list <- sort(tissue_number_list)
tissue_number_list <- rep(tissue_number_list,each=15)
                            
                            
## gene_number
files_input = list.files(path = '/gpfs/ysm/scratch60/cpsc424/wd262/sc_immune/simulation/gene_number',pattern = "*", full.names = TRUE)
gene_number_list <- sapply(files_input, basename)
files_input <- as.vector(sapply(files_input, function(file) paste0(file,'/',c(101:115),'/')))
files_W_empirical = as.vector(sapply(files_input, function(file) paste0(file,'W_empirical.txt')))
files_W_transig = as.vector(sapply(files_input, function(file) paste0(file,'W_transig.txt')))
files_Y0 = as.vector(sapply(files_input, function(file) paste0(file,'Y0.txt')))  
                                    
names(gene_number_list) <- NULL
gene_number_list <- as.numeric(gene_number_list)
#tissue_number_list <- sort(tissue_number_list)
gene_number_list <- rep(gene_number_list,each=15) 
                            
                            
## unbalanced
files_rst = "/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/simulation_rst/unbalanced/rst.T=5.D=500.K=8.corrupt=0.25.tauW=1.tauXdBeta=1.rds"
files_rst_frac = "/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/simulation_rst/unbalanced/rst_frac.T=5.D=500.K=8.corrupt=0.25.tauW=1.tauXdBeta=1.rds"
files_input = paste0('/gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/unbalanced/')
files_W_empirical = as.vector(sapply(files_input, function(file) paste0(file,'W_empirical.txt')))
files_W_transig = as.vector(sapply(files_input, function(file) paste0(file,'W_transig.txt')))
files_Y0 = as.vector(sapply(files_input, function(file) paste0(file,'Y0.txt'))) 
                            
                            
## ggplot p-val plotting
library(ggpubr)
library(farver)
compare_means(Neutrophil ~ Patient_Group, data = metadata_cnt)
my_comparisons <- list( c('ICU','Non-ICU'), c('Non-ICU','Recovery') , c('ICU','Recovery'))
p <- ggboxplot(metadata_cnt, x = "Patient_Group", y = "Neutrophil",
          color = "Patient_Group", palette = "jco",
          add = "jitter")
pdf("boxplot_Neutrophil.pdf")
p + stat_compare_means(comparisons = my_comparisons)
#p
dev.off()                            