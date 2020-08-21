same_input$X[[1]] <- readRDS('data/Pseudo_Bulk/PBMC_ImmuneCell_updated.rds')
same_input$X[[2]] <- readRDS('data/Pseudo_Bulk/Lung_ImmuneCell_updated.rds')
saveRDS(same_input,'data/Pseudo_Bulk/SAME_Input.rds')



W_T <- w_est
pi_pre <- pi_ber_est[j]
v <- v_est[[Lambda[i+1]]]
tau_w <- tau_w_est[j]


tau_x_est[j+1] <- update_tau_x(X, 
                            w_est,
                            alpha_x = alpha_prior_x, beta_x = beta_prior_x,
                            C0, c_k, K, T, D
                            )                      