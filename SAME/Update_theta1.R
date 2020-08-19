# Update MMAP parameters, theta1
set.seed(2020)
## input conventions:
## theta1=(z, w); z and w are lists of matrices
## theta2=(tau_e, alpha, gamma, v, pi, tau_x, tau_w, tau_v); each variable is a x by Gamma_i matrix
## Global params: Y0, X0, W_tilde; t0, N, D, K, C(k)=counts of celltype k, Cl=celltype lables

update_z <- function(tau_e, alpha, gamma, v, pi, tau_x, tau_w, tau_v, z, w, t0=1){
  Gammai = ncol(tau_e)
  t0 = 1
  w_t0 = w[[t0]]
  ## here below z is updated by row
  z_new <- z
  w_pseudo <- sapply(1:Gammai, function(j) alpha[1,j] * W_tilde + (1-alpha[1,j]) * w_t0) ## length of Gammai list
  tau_z <- rowSums(sapply(1:Gammai, function(j) colSums(tau_e[,j] * w_pseudo[[j]]^2) + 1))
  for(k in 1:K){
    tau_zk <- 1/tau_z[k]
    res <- lapply(w_pseudo, function(w) Y0 - w %*% z_new + w[,k] %*% z_new[k,])
    mu_zk <- rowSums(sapply(1:Gammai, function(j) colSums(tau_e[,j] * w_pseudo[[j]][,k] * res[[j]])))/tau_zk
    z_new[k,] <- sapply(mu_zk, function(mu_zkn) rnorm(mu_zkn, 1/tau_zk, 1))
  }
  return(z_new)
}

update_w <- function(tau_e, alpha, gamma, v, pi, tau_x, tau_w, tau_v, z, w, t0=1){
  Gammai = ncol(tau_e)
  for(t in 1:T){
    wt_new <- w[[t]]
    Xt <- X[[t]]
    ## here below w is updated by each element w_t0[d,k]
    if(t == t0){
      tau_wt <- sapply(1:K, function(k) rowSums(sapply(1:Gammai,
        function(j) (1 - alpha[1,j])^2 * tau_e[,j] * rowSums(z[k,]^2) + C[k] * tau_x[,j] + tau_w[1,j])))  ## D by K matrix
      for(k in 1:K){
        tau_wtk <- tau_wt[,k] ## vector of length D
        for(d in 1:D){
          res_d <- sapply(1:Gammai, function(j) sum(Y0[d,] - (alpha[1,j] * W_tilde[d,] + (1-alpha[1,j]) * wt_new[d,]) %*% z) ## vector of length Gammai
          mu_wt_dk <- sum(sapply(1:Gammai, function(j){
            tau_e[d,j] * (1 - alpha[1,j]) * sum(z[k,] * (res_d[j] + (1 - alpha[1,j]) * z[k,])) +
            tau_x[d,j] * sum(Xt[d, Cl == k]) + gamma[[j]][d,k] * v[[j]][d,k] * tau_w[1,j]))/tau_wtk[d]
          }
          wt_new[d,k] <- rnorm(mu_wt_dk, tau_wtk[d], 1)
        }
      }
    } else{
      tau_wt <- rowSums( C[k] * tau_x + t(matrix(rep(tau_w, nrow(tau_x)), nrow=Gammai))  ## d by 1
      for(k in 1:K){
        mu_wtk <- rowSums(sapply(1:Gammai, function(j) tau_x[,j] * rowSums(Xt[, Cl == k]) + gamma[[j]][,k] * v[[j]][,k] * tau_w[1,j]))/tau_wt
        wt_new[,k] <- sapply(1:D, function(d) rnorm(mu_wtk[d], tau_wt[d], 1))
      }
    }
  }
  return(wt_new)
}
