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
    ## Input:
    ## Y0: bulk sample, D*N
    ##
    ## Output:
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
                                w_t0, z, tau_alpha, tau_e
                                ){
    para1 <- sum(tau_e %*% (((W_tilde - w)%*%z)^2))
    para2 <- (1/para1) * sum(tau_e %*% ((Y0 - w_t0%*%z)*(W_tilde - w_t0)%*%z))
    alpha_unif_new <- rnorm(1, para2, para1)
    return(alpha_unif_new)
}

#update gamma
update_gamma <- function(W_tilde, 
                            W_T, pi_est, v
                            D, K
                            ){
    gamma_new <- matrix(0, nrow = D, ncol = K)
    para <- log(pi_est/(1-pi_est)) - (tau_w/2)*Reduce("+", lapply(W_T, function(x)(x-v)^2)) + (tau_w/2)*Reduce("+", lapply(W_T, function(x)x^2))
    for (d in 1:D){
        for (k in 1:K)P{
            p_temp <- 1/(1+exp(-para[d,k]))
            gamma_new[d,k] = rbinom(n=1, size = 1, prob = p_temp)
        }
    }
    return(gamma_new)
}