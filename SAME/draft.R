same_input$X[[2]] <- readRDS('data/Pseudo_Bulk/PBMC_ImmuneCell_updated.rds')
same_input$X[[1]] <- readRDS('data/Pseudo_Bulk/Lung_ImmuneCell_updated.rds')
saveRDS(same_input,'data/Pseudo_Bulk/SAME_Input.rds')



tau_e = tau_e_same
alpha = alpha_same
z = z_est
w = w_est

z_est <- update_z(Y0, X0, W_tilde,
                      tau_e_same, alpha_same, z_est, 
                      w_est, t0=1
                      )