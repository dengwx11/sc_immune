# Update nuisance parameters, theta2
# output is matrix
set.seed(2020)

update_tau_e <- function(Y0, W_tilde, # input
                            alpha, w, z, # est
                            alpha_e = alpha_prior_e, beta_e = beta_prior_e, # noninformative prior
                            N, D # other para
                            ){                                
    para1 <- rep(N/2 + alpha_e,D)
    para2 <- 0.5 * apply( (Y0 - (alpha*W_tilde+(1-alpha)*w)%*%z)^2, 1, sum) + beta_e
    tau_e_new <- matrix(0,nrow = 1, ncol = D)
    for(d in 1:D){
        tau_e_new[,d] <- rgamma(1, para1[d],para2[d])
    }
    return(tau_e_new)
}