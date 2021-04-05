set.seed(2021)
library(nnls)
library(ggplot2)
library(plotROC)
library(Seurat)
library(gridExtra)
options(repr.plot.width=14, repr.plot.height=10)

files_rst = list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/simulation_rst/corrupt',pattern = "rst.T=5*", full.names = TRUE)
files_rst_frac = list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/simulation_rst/corrupt',pattern = "rst_frac.*", full.names = TRUE)
files_same_input = paste0('/gpfs/loomis/scratch60/zhao/wd262/sc_immune/simulation/corrupt/',c(1:100),'/raw_same_input.rds')


get_summary.info <- function(files_rst,files_rst_frac,files_same_input){
    summary.info <- list()

    summary.info$cor.same <- c()
    summary.info$cor.music <- c()
    summary.info$cor.nnls <- c()

    summary.info$z_true <- list()
    summary.info$z_same <- list()
    summary.info$z_music <- list()
    summary.info$z_nnls <- list()

    summary.info$W_tilde <- list()
    summary.info$est_gamma <- list()
    summary.info$true_gamma <- list()

    summary.info$true_vg <- list()
    summary.info$est_vg <- list()

    summary.info$true_w <- list() 
    summary.info$est_w <- list() 

    mcmc_samples_theta1 =100

    for(i in seq(length(files_rst))){
    #for(i in 1:2){
        rst <- readRDS(files_rst[i])
        rst_frac <- readRDS(files_rst_frac[i])
        same_input <- readRDS(files_same_input[i])

        true_z = t(same_input$true_Z)
        true_w =  same_input$true_w$w ## true w for T tissues
        true_v = same_input$true_w$v 
        true_gamma = same_input$true_w$gamma
        W_tilde = same_input$W_tilde 


        summary.info$z_tranSig[[i]] <- rst_frac$z_est_tranSig
        summary.info$z_music[[i]] <- rst_frac$z_est_music
        summary.info$z_empirical[[i]] <- rst_frac$z_est_empirical
        summary.info$z_true[[i]] <- true_z

        summary.info$cor.same <- append(summary.info$cor.same, cor(as.vector(rst_frac$z_est_tranSig),as.vector(true_z)))
        summary.info$cor.music <- append(summary.info$cor.music, cor(as.vector(rst_frac$z_est_music),as.vector(true_z)))
        summary.info$cor.empirical <- append(summary.info$cor.empirical, cor(as.vector(rst_frac$z_est_empirical),as.vector(true_z)))


        summary.info$W_tilde[[i]] <- W_tilde
        summary.info$est_gamma[[i]] <- rst$averg_gamma
        summary.info$true_gamma[[i]] <- true_gamma

        summary.info$true_vg[[i]] <- true_v * true_gamma 
        summary.info$est_vg[[i]] <-  rst$vg

        summary.info$true_w[[i]] <- true_w 
        summary.info$est_w[[i]] <-  rst$theta1$w[[1]]

    }
    return(summary.info)
}

summary.info <- get_summary.info(files_rst,files_rst_frac,files_same_input)
saveRDS(summary.info,"/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/simulation_rst/summary.info/summary.info.corrupt.rds")

## correlation plot
corrupt_pi <- c(1:100)/100
dat <- data.frame(cor = c(summary.info$cor.same,summary.info$cor.music,summary.info$cor.empirical),
                 corrupt_pi = rep(c(1:100)/100,3),
                 method = c(rep("SAME",100),rep("MuSiC",100),rep("empirical",100)))

ggplot(dat, aes(corrupt_pi,cor)) + 
geom_point(aes(colour = method), shape = 4) +
labs(y = 'Correlation', x = "Corruption Rate", title = 'Correlation with True Cell Proportions') + 
geom_line(aes(color = method, linetype = method), size=2) +
# ylim(0.8, 1) +
theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))

## ROC plot
tbl <- lapply(c(1:7),function(idx) 
    data.frame('Estmated_Gamma' = as.vector(summary.info$est_gamma[[10*(idx-1)+1]]),
               "True_Gamma" = as.vector(summary.info$true_gamma[[10*(idx-1)+1]]),
               "corrupt" = rep((idx-1)*0.1,4000)))
tbl <- Reduce(rbind, tbl)
              
ggplot(tbl, aes(d=True_Gamma,m=Estmated_Gamma,color=factor(corrupt)))+geom_roc(n.cuts = 0)     
              
              
##scatterplot    
i=50
scatter.df <- data.frame(tranSig = as.vector(summary.info$z_tranSig[[i]]),
                         music = as.vector(summary.info$z_music[[i]]),
                         empirical = as.vector(summary.info$z_empirical[[i]]),
                         true = as.vector(summary.info$z_true[[i]]),
                         celltype = rep(paste0('celltype',c(1:8)),each = 200))
p <- list()
p[[1]]<-ggplot(scatter.df, aes(x=true, y = tranSig, color = celltype))+geom_point()+ coord_fixed(ratio = 1)+ xlim(0,1)+ ylim(0,1) 
p[[2]]<-ggplot(scatter.df, aes(x=true, y = music, color = celltype))+geom_point()+ coord_fixed(ratio = 1)+ xlim(0,1)+ ylim(0,1)+theme(legend.position="none")
p[[3]]<-ggplot(scatter.df, aes(x=true, y = empirical, color = celltype))+geom_point()+ coord_fixed(ratio = 1)+ xlim(0,1)+ ylim(0,1)+theme(legend.position="none")
legend <- get_legend(p[[1]])
p[[1]] <- p[[1]] + theme(legend.position="none")


pll <- grid.arrange(p[[1]],p[[2]],p[[3]],legend,ncol=4)              