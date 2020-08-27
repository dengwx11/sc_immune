same_input$X[[2]] <- readRDS('data/Pseudo_Bulk/PBMC_ImmuneCell_updated.rds')
same_input$X[[1]] <- readRDS('data/Pseudo_Bulk/Lung_ImmuneCell_updated.rds')
saveRDS(same_input,'data/Pseudo_Bulk/SAME_Input.rds')



tau_e <- tau_e_same
alpha <- alpha_same
z <- z_est
w <- w_est
gamma <-gamma_same
v <- v_same
tau_x <- tau_x_same
tau_w <- tau_w_same

z_est <- update_z(Y0, X0, W_tilde,
                      tau_e, alpha, z, w, t0=1)

w_p <-w_pseudo[[1]]