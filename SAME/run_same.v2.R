set.seed(2021)
source("SAME/Update_theta1.v2.R")
source("SAME/Update_theta2.v2.R")
source("SAME/clean_format.R")


#### demo
#setwd("~/zhao-data/sc_immune/sc_immune/")
#source('Batch_Correction/get_tissue_specific_input.R')
#source('SAME/run_same.v2.R')
#source("SAME/clean_format.R")

#taget_tissue <- 'PBMC'
#files = list.files(path = '/gpfs/ysm/home/bl666/HCL/Pseudo_Bulk',pattern = "*_updated.rds", full.names = TRUE)
#tissue_list <- c("BM","CB","Kidney","Liver","Lung","PBMC")
#celltype_used_list <- c("B", "Neutrophil", "NK cell", "T")
#YSG <- readRDS("/gpfs/ysm/pi/zhao-data/wd262/sc_immune/sc_immune/data/NSCLC/sg.list.rds")
#output_path = "/gpfs/ysm/pi/zhao-data/wd262/sc_immune/write/pipeline_on_HCL"

#rst <- run_SAME(target_tissue,tissue_list,celltype_used_list,files,YSG,empirical_pi=0.3,liger.turnon=FALSE,output_path)


run_SAME <- function(target_tissue,tissue_list,celltype_used_list,files,YSG,empirical_pi,liger.turnon=TRUE,output_path){
    
    ans <- get_tissue_specific_input(target_tissue,tissue_list,celltype_used_list,files,YSG,output_path,liger.turnon=liger.turnon)
    
    YSG <- ans$YSG
    w_tissue_indicator <- ans$w_tissue_indicator
    seur.TPM_list <- ans$seur.TPM_list
    seur_list <- ans$X
    liger <- ans$liger
    
    celltype_lists <- lapply(seq(length(seur_list)),function(i) seur_list[[i]]$Celltype_used)
    used_cell <- get_used_cell(celltype_lists,celltype_used_list)
    ans <- keep_used_cell(seur_list,used_cell, YSG)
    X_mat <- ans$X_mat  
    celltype_lists <- ans$celltype_lists    
        
        
    c_k<-get_c_k(celltype_lists,celltype_used_list)
    
    
                             
                             
    Cl <- get_Cl(celltype_lists,celltype_used_list)$Cl 
                             
    
                             
    K = length(celltype_used_list)
    T = length(tissue_list)
    D = length(YSG)                         
    rst <- SAME(X_mat,w_tissue_indicator,Cl,c_k,K,T,D,
                 empirical_pi,
                 mcmc_samples_theta1, Lambda, YSG)
    rst$w_tissue_indicator <- w_tissue_indicator                         
    return(rst)                         
}

SAME <- function(X_mat,w_tissue_indicator,Cl,c_k,K,T,D,
                 empirical_pi,
                 mcmc_samples_theta1, Lambda, YSG)
{
    
    ##### Input ######
    ## X_mat: list of expr matrix
    ## w_tissue_indicator: list of tissue-celltype-specific gene matrix T(tissue)*D(gene)*K(celltype)
    ## empiricla_pi
    
    
  

  
  
  ## Initialization (lower case)
  # Noninformative Prior Parameters
  
  tau_v = 0.6
  alpha_prior_e = 10^(-6)
  beta_prior_e = 10^(-6)
  alpha_prior_x = 10^(-6)
  beta_prior_x = 10^(-6)
  alpha_prior_w = 10^(-6)
  beta_prior_w = 10^(-6)
  alpha_prior_pi = 10^(-6)
  beta_prior_pi = 10^(-6)
  
  
  
  # Estimations
  # numerber or vector is stored as mcmc sequence
  # matrix is stored as the newest one 
  # theta2
  mcmc_samples_theta2 = sum(Lambda)

  gamma_est <- list()
  for(i in 1:Lambda[length(Lambda)]){
    gamma_est[[i]] <- matrix(rbinom(D*K,1,0.5),nrow = D, ncol = K)
    #gamma_est[[i]] <- true_gamma
  }
  v_est <- list()
  for(i in 1:Lambda[length(Lambda)]){
    v_est[[i]] <- abs(matrix(rnorm(D*K,mean = 2,sd = 1),nrow = D, ncol = K))
    # v_est[[i]] <- true_v
  }
  
  pi_ber_est <- matrix(empirical_pi, nrow =  mcmc_samples_theta2+1, ncol = 1)
  pi_ber_est[1] <- empirical_pi
  tau_x_est <- matrix(1000, nrow =  mcmc_samples_theta2+1, ncol = D)
  tau_x_est[1,] <- 1
  #tau_x_est <- matrix(raw_X[[1]]$tau_xd, nrow =  1, ncol = D)[rep(1,mcmc_samples_theta2+1),]
  tau_w_est <- matrix(0, nrow =  mcmc_samples_theta2+1, ncol = 1)
  tau_w_est[1] <- 10
  
  # theta1

  w_est <- list()
  for(i in 1:T) w_est[[i]] <- matrix(rbinom(D*K,1,0.5),nrow = D, ncol = K) # first item = t0 tissue 

  #w_est = same_input$true_w$w
  
  
  ## SAME
  ptm <- proc.time()
  for(i in 1:mcmc_samples_theta1){
    ## Step 1: Update theta2 for Lambda[i] times
    print(Lambda[1:(i+1)])
    start_idx = sum(Lambda[1:i])+1
    end_idx = sum(Lambda[1:(i+1)])        
    Lambdai = Lambda[i+1]
    print(c(start_idx,end_idx, Lambdai))
    
    k=1
    j = start_idx ## for debugging
    
    for(j in start_idx:end_idx){
      

      
    
      if( i  == 1) {
        last_v = 1
      }else if( k == 1) 
      {
        last_v = Lambda[i+1]
      }else {last_v = k-1}
      gamma_est[[k]] <- update_gamma(w_est, pi_ber_est[j], v_est[[ last_v ]], tau_w_est[j],w_tissue_indicator)

      
      v_est[[k]]<-update_v(tau_w_est[j], w_est, gamma_est[[k]],w_tissue_indicator,
                           tau_v=tau_v
      )

        
      # pi_ber_est[j+1] <- update_pi_est(gamma_est[[k]], 
      #                                  alpha_pi = alpha_prior_pi, beta_pi = beta_prior_pi
      # )
      
      tau_x_est[j+1,] <- update_tau_x(X_mat, w_tissue_indicator,
                      w_est, Cl, c_k,
                      alpha_x = alpha_prior_x, beta_x = beta_prior_x
                      )
                       
      tau_w_est[j+1] <- update_tau_w(w_est, v_est[[k]], gamma_est[[k]],w_tissue_indicator,
                                     alpha_w = alpha_prior_w, beta_w = beta_prior_w
      )   
                
      k=k+1
    }
    
    print(table(gamma_est[[k-1]]))
    
    
    ## Step 2: Update theta1 for one time
    gamma_same <- gamma_est[1:(k-1)]
    v_same <- v_est[1:(k-1)]
    tau_x_same <- t(matrix(tau_x_est[(start_idx+1):(end_idx+1),],nrow = Lambdai, ncol = D))
    tau_w_same <- matrix(tau_w_est[(start_idx+1):(end_idx+1),], nrow = 1, ncol = Lambdai)

    
    w_est <- update_w(X_mat, 
                      gamma_same, v_same, tau_x_same, tau_w_same, w_est, Cl
    )
    
  }
  print(proc.time() - ptm)
  
  rst <- list()
  rst$theta1 <- list()
  rst$theta2 <- list()
  

  rst$theta1$w <- w_est
  

  rst$theta2$gamma <- gamma_est
  rst$theta2$v <- v_est
  rst$theta2$pi_ber <- pi_ber_est
  rst$theta2$tau_x <- tau_x_est
  rst$theta2$tau_w <- tau_w_est
    
  ans <- get_est_vg(rst)
  rst$vg <-  ans$vg
  rst$averg_gamma <- ans$averg_gamma  
  
  return(rst)
  
}
