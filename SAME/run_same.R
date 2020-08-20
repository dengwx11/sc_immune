set.seed(2020)
source("Update_theta1.R")
source("Update_theta2.R")

SAME <- function(Y0, X, W_tilde,
                 mcmc_samples_theta1, Lambda)
{


    ## Initialization (lower case)
    # Noninformative Prior Parameters

    tau_v = 0.01
    alpha_prior_e = 10^(-6)
    beta_prior_e = 10^(-6)
    alpha_prior_x = 10^(-6)
    beta_prior_x = 10^(-6)
    alpha_prior_w = 10^(-6)
    beta_prior_w = 10^(-6)
    alpha_prior_pi = 10^(-6)
    beta_prior_pi = 10^(-6)


    D = nrow(Y0) # number of genes
    K = ncol(W_tilde) # number of cell types
    T = length(X)  # number of tissues
    N = ncol(Y0) # number of bulk samples
    c_k = matrix(200,nrow = K, ncol = T) # sequence of cell counts for each cell type and each tissue
    C0 = sum(c_k) # total number of single cell counts


    # Estimations
    # numerber or vector is stored as mcmc sequence
    # matrix is stored as the newest one 
    # theta2
    mcmc_samples_theta2 = sum(Lambda)
    tau_e_est <- matrix(0, nrow =  mcmc_samples_theta2, ncol = D)
    tau_e_est[1] <- 1
    alpha_unif_est <- matrix(0, nrow =  mcmc_samples_theta2, ncol = 1)
    alpha_unif_est[1] <- 0.5
    gamma_est <- list()
    for(i in 1:Lambda[length(Lambda)]){
        gamma_est[[i]] <- matrix(rbinom(D*K,1,0.5),nrow = D, ncol = K)
    }
    v_est <- list()
    for(i in 1:Lambda[length(Lambda)]){
        v_est[[i]] <- matrix(rnorm(D*K,mean = 2,sd = 1),nrow = D, ncol = K)
    }
    pi_ber_est <- matrix(0, nrow =  mcmc_samples_theta2, ncol = 1)
    pi_ber_est[1] <- 0.5
    tau_x_est <- matrix(0, nrow =  mcmc_samples_theta2, ncol = D)
    tau_x_est[1] <- 1
    tau_w_est <- matrix(0, nrow =  mcmc_samples_theta2, ncol = 1)
    tau_w_est[1] <- 1

    # theta1
    z_est <- matrix(rnorm(K*N,mean = 2,sd = 1),nrow = K, ncol = N)
    w_est <- list()
    for(i in 1:T) w_est[[i]] <- matrix(rbinom(D*K,1,0.5),nrow = D, ncol = K) # first item = t0 tissue 


    ## SAME

    for(i in 1:mcmc_samples_theta1){
        ## Step 1: Update theta2 for Lambda[i] times
        start_idx = sum(Lambda[1:i])
        end_idx = sum(Lambda[1:(i+1)])
        k=1
        for(j in start_idx:end_idx){
            tau_e_est[j+1] <- update_tau_e(Y0, W_tilde, # input
                            alpha[j], w_t0[[1]], z_est, # est
                            alpha_e = alpha_prior_e, beta_e = beta_prior_e, # noninformative prior
                            N, D # other para
                            )
            alpha_unif_est[j+1] <- update_alpha_unif(Y0, W_tilde,
                                w_t0[[1]], z_est, tau_e_est[j+1]
                                )
            gamma_est[[k]] <- update_gamma(W_tilde, 
                            w_est, pi_ber_est[j], v_est[[Lambda[i+1]]],
                            D, K
                            )
            v_est[[k]]<-update_v(tau_w, w_est, gamma_est[[k]],
                        tau_v=tau_v,
                        D, K, T
                        )
            pi_ber_est[j+1] <- update_pi_est(gamma_est[[k]], 
                            alpha_pi = alpha_prior_pi, beta_pi = beta_prior_pi,
                            D, K
                            )
            tau_x_est[j+1] <- update_tau_x(X, 
                            w_est,
                            alpha_x = alpha_prior_x, beta_x = beta_prior_x,
                            C0, c_k, K, T, D
                            )
            tau_w_est[j+1] <- update_tau_w(w_est, v_est[[k]], gamma_est[[k]],
                            alpha_w = alpha_prior_w, beta_w = beta_prior_w,
                            T, D, K
                            )                
            k=k+1
        }

        ## Step 2: Update theta1 for one time
        z_est <- update_z(theta)
        w_est <- update_w(theta)
    }

    rst <- list()
    rst$theta1 <- list()
    rst$theta2 <- list()

    rst$theta1$z <- z_est
    rst$theta1$w <- w_est

    rst$theta2$tau_e <- tau_e_est
    rst$theta2$alpha_unif <- alpha_unif_est
    rst$theta2$gamma <- gamma_est
    rst$theta2$v <- v_est
    rst$theta2$pi <- pi_ber_est
    rst$theta2$tau_x <- tau_x_est
    rst$theta2$tau_w <- tau_w_est

    return(rst)

}