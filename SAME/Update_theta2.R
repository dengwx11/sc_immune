# Update nuisance parameters, theta2
set.seed(2020)

update_tau_e <- function(Y0, W_tilde, alpha, w, z, alpha_e = alpha_prior_e, beta_e = beta_prior_e, N, D, Lambda){
    para1 <- rep(N/2 + alpha_e,D)
    para2 <- 0.5 * apply( Y0 - (alpha*w_tilde+(1-alpha)*w)%*%z, 1, sum) + beta_e
    tau_e_new <- matrix(0,nrow = Lambda, ncol = D)
    for(d in 1:D){
        tau_e_new[,d] <- rgamma(Lambda, para1[d],para2[d])
    }
    return(tau_e_new)
}