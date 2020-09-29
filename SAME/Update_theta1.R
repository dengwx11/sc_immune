# Update MMAP parameters, theta1
set.seed(2020)
## input conventions:
## theta1=(z, w); z and w are lists of matrices
## theta2=(tau_e, alpha, gamma, v, pi, tau_x, tau_w, tau_v); each variable is a x by Gamma_i matrix
## Global params: Y0, X0, W_tilde; t0, N, D, K, C(k)=counts of celltype k, Cl=celltype lables

update_z <- function(Y0, X0, W_tilde,
                      tau_e, alpha, z, w, t0=1){
 ####
 ## vector -> matrix ( D * Lambdai )
 ## matrix -> list (length = Lambdai )
 ####                        
  Lambdai = ncol(tau_e)
  t0 = 1
  w_t0 = w[[t0]]
  ## here below z is updated by row
  z_new <- z
  w_pseudo <- lapply(1:Lambdai, function(j) alpha[1,j] * W_tilde + (1-alpha[1,j]) * w_t0) ## length of Lambdai list

  tau_z <- sapply(1:Lambdai, function(j) colSums(tau_e[,j] * w_pseudo[[j]]^2) + 1 )
  tau_z <- matrix(tau_z, nrow = Lambdai, ncol = K)
  
  for(k in 1:K){
    tau_zk <- tau_z[k]
    res <- lapply(w_pseudo, function(w_p) Y0 - w_p %*% z_new )
    mu_zk <- sapply(1:Lambdai, function(j) tau_e[,j] * matrix(w_pseudo[[j]][,k],nrow = 1) %*% res[[j]])/tau_zk
    mu_zk <- matrix(mu_zk, nrow = Lambdai, ncol = N)
    z_new[k,] <- sapply(mu_zk, function(mu_zkn) rnorm(1, mu_zkn, 1/tau_zk))
  }

  return(z_new)
  
}

update_w <- function(Y0, X0, W_tilde,
                      tau_e, alpha, gamma, v, tau_x, tau_w, z, w, t0=1
                      ){
  Lambdai = ncol(tau_e)
  end_idx = 0
  for(t in 1:T){
    wt_new <- w[[t]]
    Xt <- X[[t]]
    ## here below w is updated by each element w_t0[d,k]
    C = as.vector(c_k[,t])
    start_idx <- end_idx + 1
    end_idx <- start_idx + sum(c_k[,t])-1
    if(t == t0){

              tau_wt <- sapply(1:K, function(k) sapply(1:Lambdai,
        function(j) (1 - alpha[1,j])^2 * tau_e[,j] * z[k,]^2 + C[k] * tau_x[,j] + tau_w[1,j]))  ## D by K matrix
        tau_wt <- matrix(tau_wt, nrow = Lambdai, ncol = K)
        tau_wt <- apply(tau_wt, 1, sum)

      for(k in 1:K){
        tau_wtk <- tau_wt[,k] ## vector of length D
        for(d in 1:D){
          res_d <- sapply(1:Lambdai, function(j) sum(Y0[d,] - (alpha[1,j] * W_tilde[d,] + (1-alpha[1,j]) * wt_new[d,]) %*% z))## vector of length Lambdai
          mu_wt_dk <- sum(sapply(1:Lambdai, function(j){
            (tau_e[1,j] * (1 - alpha[1,j]) * sum(z[k,] * (res_d[j] + (1 - alpha[1,j]) * z[k,])) +
            tau_x[d,j] * sum(Xt[d, Cl[start_idx:end_idx] == k]) + gamma[[j]][d,k] * v[[j]][d,k] * tau_w[1,j])/tau_wtk[d]
          }))
          wt_new[d,k] <- rnorm(mu_wt_dk, tau_wtk[d], 1)
        }
      }
    } else{
      tau_wt <- rowSums( C[k] * tau_x + t(matrix(rep(tau_w, nrow(tau_x)), nrow=Lambdai)))  ## d by 1
      for(k in 1:K){
        mu_wtk <- rowSums(sapply(1:Lambdai, function(j) tau_x[,j] * rowSums(Xt[, Cl[start_idx:end_idx] == k]) + gamma[[j]][,k] * v[[j]][,k] * tau_w[1,j]))/tau_wt
        wt_new[,k] <- sapply(1:D, function(d) rnorm(mu_wtk[d], tau_wt[d], 1))
      }
    }
  }
  return(wt_new)
}


