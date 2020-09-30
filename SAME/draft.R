same_input$X[[2]] <- readRDS('data/Pseudo_Bulk/PBMC_ImmuneCell_updated.rds')
same_input$X[[1]] <- readRDS('data/Pseudo_Bulk/Lung_ImmuneCell_updated.rds')
saveRDS(same_input,'data/Pseudo_Bulk/SAME_Input.rds')

# debug tau_e

alpha = alpha_unif_est[j]
w_t0 = w_est[[1]]
z = z_est

# debug gamma
W_T <- w_est
pi_pre <- pi_ber_est[j]
v <- v_est[[ last_v ]]
tau_w <- tau_w_est[j]

# debug v
tau_w <- tau_w_est[j]
W_T <- w_est
gamma <- gamma_est[[k]]


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


tau_e = 100
W = (1 - alpha)*(true_w_sc[[1]]) + alpha*(W_tilde)
    e <- matrix(rnorm(D*N, 0, sd = (1/sqrt(tau_e))), nrow = D, ncol = N)
    e[e < 0] = 0.001
    Y = W %*% Z + e                      