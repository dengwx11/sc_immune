same_input$X[[2]] <- readRDS('data/Pseudo_Bulk/PBMC_ImmuneCell_updated.rds')
same_input$X[[1]] <- readRDS('data/Pseudo_Bulk/Lung_ImmuneCell_updated.rds')
saveRDS(same_input,'data/Pseudo_Bulk/SAME_Input.rds')

# debug Z
tau_e = tau_e_same
alpha = alpha_same
z = z_est
w = w_est







# debug W
tau_e = tau_e_same
alpha = alpha_same
gamma = gamma_same
v = v_same
tau_x =  tau_x_same
tau_w = tau_w_same
z = z_est
w = w_est


w_est <- update_w(Y0, X0, W_tilde,
                      tau_e_same, alpha_same, gamma_same, v_same, tau_x_same, tau_w_same, z_est, w_est, t0=1
                      )