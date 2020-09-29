same_input$X[[2]] <- readRDS('data/Pseudo_Bulk/PBMC_ImmuneCell_updated.rds')
same_input$X[[1]] <- readRDS('data/Pseudo_Bulk/Lung_ImmuneCell_updated.rds')
saveRDS(same_input,'data/Pseudo_Bulk/SAME_Input.rds')



W_tilde <- W_tilde
W_T <- w_est
pi_pre <- pi_ber_est[j]
v <- v_est[[ last_v ]]
tau_w <- tau_w_est[j]

update_gamma(W_tilde, W_T, pi_pre, v, tau_w)