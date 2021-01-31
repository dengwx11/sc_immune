set.seed(2020)
source("SAME/Update_theta1.R")
source("SAME/Update_theta2.R")
source("SAME/clean_format.R")





SAME <- function(Y0, X, W_tilde,
                 mcmc_samples_theta1, Lambda, c_k, YSG, alpha = 0.5)
{
  
  Cl <- get_Cl(X)
  celltype_list <- Cl$levels
  Cl = Cl$Cl
  
  X_mat <- get_X_mat(X, YSG)
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
  
  
  
  # Estimations
  # numerber or vector is stored as mcmc sequence
  # matrix is stored as the newest one 
  # theta2
  mcmc_samples_theta2 = sum(Lambda)
  tau_e_est <- matrix(1000, nrow =  mcmc_samples_theta2+1, ncol = D)
  tau_e_est[1] <- 1000
  alpha_unif_est <- matrix(alpha, nrow =  mcmc_samples_theta2+1, ncol = 1)
  #alpha_unif_est[1] <- 0.5
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
  
  pi_ber_est <- matrix(0.3, nrow =  mcmc_samples_theta2+1, ncol = 1)
  pi_ber_est[1] <- 0.3
  tau_x_est <- matrix(1000, nrow =  mcmc_samples_theta2+1, ncol = D)
  tau_x_est[1,] <- 1
  #tau_x_est <- matrix(raw_X[[1]]$tau_xd, nrow =  1, ncol = D)[rep(1,mcmc_samples_theta2+1),]
  tau_w_est <- matrix(0, nrow =  mcmc_samples_theta2+1, ncol = 1)
  tau_w_est[1] <- 10
  
  # theta1
  z_est <- abs(matrix(rnorm(K*N,mean = 0.8,sd = 1),nrow = K, ncol = N))
  w_est <- list()
  for(i in 1:T) w_est[[i]] <- matrix(rbinom(D*K,1,0.5),nrow = D, ncol = K) # first item = t0 tissue 
  #z_est = true_z
  #w_est = same_input$true_w$w
  
  
  ## SAME
  ptm <- proc.time()
  for(i in 1:mcmc_samples_theta1){
    ## Step 1: Update theta2 for Lambda[i] times
    start_idx = sum(Lambda[1:i])+1
    end_idx = sum(Lambda[1:(i+1)])        
    Lambdai = Lambda[i+1]
    print(c(start_idx,end_idx, Lambdai))
    
    k=1
    j = start_idx ## for debugging
    
    for(j in start_idx:end_idx){
      tau_e_est[j+1,] <- update_tau_e(Y0, W_tilde, # input
                                      alpha_unif_est[j], w_est[[1]], z_est, # est
                                      alpha_e = alpha_prior_e, beta_e = beta_prior_e # noninformative prior
      )
      #alpha_unif_est[j+1,] <- update_alpha_unif(Y0, W_tilde,
      #                    w_est[[1]], z_est, tau_e_est[j+1,]
      #                    )
      if( i  == 1) {
        last_v = 1
      }else if( k == 1) 
      {
        last_v = Lambda[i+1]
      }else {last_v = k-1}
      gamma_est[[k]] <- update_gamma(W_tilde, 
                                     w_est, pi_ber_est[j], v_est[[ last_v ]], tau_w_est[j]
      )                    
      
      
      v_est[[k]]<-update_v(tau_w_est[j], w_est, gamma_est[[k]],
                           tau_v=tau_v
      )
      pi_ber_est[j+1] <- update_pi_est(gamma_est[[k]], 
                                       alpha_pi = alpha_prior_pi, beta_pi = beta_prior_pi
      )
      tau_x_est[j+1,] <- update_tau_x(X, 
                      w_est,
                      alpha_x = alpha_prior_x, beta_x = beta_prior_x
                      )
      tau_w_est[j+1] <- update_tau_w(w_est, v_est[[k]], gamma_est[[k]],
                                     alpha_w = alpha_prior_w, beta_w = beta_prior_w
      )                
      k=k+1
    }
    
    print(table(gamma_est[[k-1]]))
    
    
    ## Step 2: Update theta1 for one time
    tau_e_same <- matrix(tau_e_est[(start_idx+1):(end_idx+1),], nrow = D, ncol = Lambdai)
    alpha_same <- matrix(alpha_unif_est[(start_idx+1):(end_idx+1),], nrow = 1, ncol = Lambdai)
    gamma_same <- gamma_est[1:(k-1)]
    v_same <- v_est[1:(k-1)]
    tau_x_same <- t(matrix(tau_x_est[(start_idx+1):(end_idx+1),],nrow = Lambdai, ncol = D))
    tau_w_same <- matrix(tau_w_est[(start_idx+1):(end_idx+1),], nrow = 1, ncol = Lambdai)
    
    z_est <- update_z(Y0, X_mat, W_tilde,
                      tau_e_same, alpha_same, z_est, w_est, t0=1
    )
    
    w_est <- update_w(Y0, X_mat, W_tilde,
                      tau_e_same, alpha_same, gamma_same, v_same, tau_x_same, tau_w_same, z_est, w_est, t0=1, Cl
    )
    
  }
  print(proc.time() - ptm)
  
  rst <- list()
  rst$theta1 <- list()
  rst$theta2 <- list()
  
  rst$theta1$z <- z_est
  rst$theta1$w <- w_est
  
  rst$theta2$tau_e <- tau_e_est
  rst$theta2$alpha_unif <- alpha_unif_est
  rst$theta2$gamma <- gamma_est
  rst$theta2$v <- v_est
  rst$theta2$pi_ber <- pi_ber_est
  rst$theta2$tau_x <- tau_x_est
  rst$theta2$tau_w <- tau_w_est
  
  return(rst)
  
}
