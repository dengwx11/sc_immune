# Update nuisance parameters, theta2
# output is matrix
set.seed(2020)
#update tau_e
update_tau_e <- function(Y0, W_tilde, # input
                            alpha, w_t0, z, # est
                            alpha_e = alpha_prior_e, beta_e = beta_prior_e, # noninformative prior
                            N, D # other para
                            ){
    #############
    #### Input:
    ## Y0: bulk sample, D*N
    ## W_tilde: empirical weight matrix, D*K
    ## w_t0: short for w_hatt0, D*K
    ## z: cluster fractions matrix, K*N
    #### Output:
    ## tau_e_new: 1*D
    #############                                                                    
    para1 <- rep(N/2 + alpha_e,D)
    para2 <- 0.5 * apply( (Y0 - (alpha*W_tilde+(1-alpha)*w_t0)%*%z)^2, 1, sum) + beta_e
    tau_e_new <- matrix(0,nrow = 1, ncol = D)
    for(d in 1:D){
        tau_e_new[,d] <- rgamma(1, para1[d],para2[d])
    }
    return(tau_e_new)
}

#update alpha_unif
update_alpha_unif <- function(Y0, W_tilde,
                                w_t0, z, tau_e
                                ){
    #############
    #### Input:
    ## Y0: bulk sample, D*N
    ## W_tilde: empirical weight matrix, D*K
    ## w_t0: short for w_hatt0, D*K
    ## z: cluster fractions matrix, K*N
    ## tau_e: 1*D
    #### Output:
    ## alpha_unif_new: one value falls into interval of (0, 1)
    ############# 
    para1 <- sum(tau_e %*% (((W_tilde - w_t0)%*%z)^2))
    para2 <- (1/para1) * sum(tau_e %*% ((Y0 - w_t0%*%z)*(W_tilde - w_t0)%*%z))
    alpha_unif_new <- rnorm(1, para2, sd = (1/sqrt(para1))
    return(alpha_unif_new)
}

#update gamma
update_gamma <- function(W_tilde, 
                            W_T, pi_pre, v
                            D, K
                            ){
    #############
    #### Input:
    ## W_tilde: empirical weight matrix, D*K
    ## W_T: a list of W_hat from each tissue, and the length of the list is T, the dimension of each W_hat is D*K
    ## pi_pre: estimated pi from the last step. 
    ## v: D*K
    #### Output:
    ## gamma_new: a matrix of Bernoulli variables, D*K
    #############                                 
    gamma_new <- matrix(0, nrow = D, ncol = K)
    para <- log(pi_pre/(1-pi_pre)) - (tau_w/2)*Reduce("+", lapply(W_T, function(x)(x-v)^2)) + (tau_w/2)*Reduce("+", lapply(W_T, function(x)x^2))
    for (d in 1:D){
        for (k in 1:K){
            p_temp <- 1/(1+exp(-para[d,k]))
            gamma_new[d,k] = rbinom(n=1, size = 1, prob = p_temp)
        }
    }
    return(gamma_new)
}

#update v
update_v <- function(tau_w, W_T, gamma,
                        tau_v=tau_v,
                        D, K, T
                        ){
    #############
    #### Input:
    ## tau_w: variance of W_T, all w share the same tau_w
    ## W_tilde: empirical weight matrix, D*K
    ## W_T: a list of W_hat from each tissue, and the length of the list is T, the dimension of each W_hat is D*K
    ## gamma: a matrix of Bernoulli variables, D*K 
    #### Output:
    ## v_new: D*K
    #############                              
    v_new <- matrix(0, nrow = D, ncol = K)
    para1 <- T*tau_w + tau_v
    para2 <- (1/para1)*tau_w*Reduce("+", W_T)
    for (d in 1:D){
        for (k in 1:K){
            if (gamma[d,k] == 1){
                v_new[d,k] <- rnorm(n=1, para2, sd = (1/sqrt(para1))
            }
            if (gamma[d,k] == 0){
                v_new[d,k] <- rnorm(n=1, 0, sd = (1/sqrt(tau_v)))
            }
        }
    }
    return(v_new)
}

#update pi_est
update_pi_est <- function(gamma, 
                            alpha_pi = alpha_prior_pi, beta_pi = beta_prior_pi,
                            D, K
                            ){                                                                
    para1 <- alpha_pi + sum(gamma)
    para2 <- beta_pi + D*K - sum(gamma)
    pi_new <- rbeta(n=1, para1, para2)
    return(pi_new)
}

#update tau_x
update_tau_x <- function(X, #X is a list of gene expression matrix from each tissue. For each matrix, cells should be ordered by cell types
                            W_T,
                            alpha_x = alpha_prior_x, beta_x = beta_prior_x,
                            C0, c_k, K, T, D
                            ){
    #############
    #### Input:
    ## X: a list of scRNA-seq gene expression matrix from each tissue. the dimension of X[[t]] is D*sum(c_k[,t])
    ## W_T: a list of W_hat from each tissue, and the length of the list is T, the dimension of each W_hat is D*K
    ## c_k: sequence of cell counts for each cell type and each tissue, K*T
    #### Output:
    ## tau_x_new: 1*D
    #############                                 
    tau_x_new <- matrix(0,nrow = 1, ncol = D)
    para2_list <- list()
    para2_temp <- matrix(0, nrow = 1, ncol = D)
    for (t in 1:T){
        w_t <- W_T[[t]]
        w_temp <- matrix(0, nrow = D, ncol = sum(c_k[,t]))
        for (k in 1:K){
            start_index <- sum(c_k[1:(k-1),t]) + 1 
            end_index <- sum(c_k[1:k,t)
            w_temp[,start_idx:end_idx] <- w_t[,k]
        }
        x_t <- X[[t]]
        para2_term1 <- ((x_t - w_t)^2)/2
        para2_list[[t]] <- apply(para2_term1, 1, sum)
        para2_temp = para2_temp + para2_list[[t]]
    }
    para1 <- rep(((C0/2) + alpha_x), D)
    para2 <- para2_temp + beta_x
    for(d in 1:D){
        tau_x_new[,d] <- rgamma(1, para1[d],para2[d])
    }
    return(tau_x_new)
}

#update tau_w
update_tau_w <- function (W_T, v, gamma,
                            alpha_w = alpha_prior_w, beta_w = beta_prior_w,
                            T, D, K
                            ){
    #############
    #### Input:
    ## W_T: a list of W_hat from each tissue, and the length of the list is T, the dimension of each W_hat is D*K
    ## v: D*K
    ## gamma: a matrix of Bernoulli variables, D*K 
    #### Output:
    ## tau_w_new: all the w share this same tau_w
    #############                                   
    para1 <- (T*D*K/2) + alpha_w
    para2 <- sum(Reduce("+", lapply(W_T, function(x)(x-v*gamma)^2)))/2 + beta_w
    tau_w_new <- rgamma(1, para1, para2)
    return(tau_w_new)
}