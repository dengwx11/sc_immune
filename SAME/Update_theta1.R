# Update MMAP parameters, theta1
set.seed(2020)
# input for Update_theta1.R is theta=list(theta1, theta2)
# theta2=(tau_e, alpha, gamma_kn, v_dk, pi, tau_x, tau_w, tau_v)
# theta1=(z, w) # w_t0 = w[[t0]]; t0 = 1
# w_t is short for w_hat

update_z <- function(theta){
  theta1 <- theta$theta1
  theta2 <- theta$theta2
  tau_e = theta2$tau_e ## d by Gammai matrix
  alpha = theta2$alpha ## 1 by Gammai matrix
  # gamma_kn = theta2$gamma_kn
  # v_dk = theta2$v_dk
  # pi_var = theta2$pi
  # tau_x = theta2$tau_xd
  # tau_w = theta2$tau_w
  # tau_v = theta2$tau_v
  z = theta1$z
  w = theta1$w
  Gammai = ncol(tau_e)
  t0 = 1
  w_t0 = w[[t0]]
  ## here below z is updated by row and progressively
  z_new <- z
  w_pseudo <- sapply(1:Gammai, function(j) alpha[1,j] * w_tildet0 + (1-alpha[1,j]) * w_t0 ## length of Gammai list
  tau_z <- rowSums(sapply(1:Gammai, function(j) colSums(tau_e[,j] * w_pseudo[[j]]^2) + 1))
  for(k in 1:K){
    tau_zk <- 1/tau_z[k]
    res <- lapply(w_pseudo, function(w) Y - w %*% z_new + w[,k] %*% z_new[k,])
    mu_zk <- rowSums(sapply(1:Gammai, function(j) colSums(tau_e[,j] * w_pseudo[[j]][,k] * res[[j]])))/tau_zk
    z_new[k,] <- sapply(mu_zk, function(mu_zkn) rnorm(mu_zkn, 1/tau_zk, 1))
  }
  theta1_new <- list(w=w, z=z_new)
  return(list(theta1=theta1_new, theta2=theta2))
}

update_w <- function(theta){
  theta1 <- theta$theta1
  theta2 <- theta$theta2
  tau_e = theta2$tau_e
  alpha = theta2$alpha
  gamma = theta2$gamma ## list of matrices d by k with length Gammai
  v = theta2$v
  pipi = theta2$pi
  tau_x = theta2$tau_xd ## d by Gammai matrix
  tau_w = theta2$tau_w
  tau_v = theta2$tau_v
  z = theta1$z
  w = theta1$w
  Gammai = ncol(tau_e)
  t0 = 1
  for(t in 1:T){
    if(t == t0){
      wt_new <- w[[t]]
      tau_wt <- sapply(1:K, function(k) rowSums(sapply(1:Gammai,
        function(j) (1 - alpha[1,j])^2 * tau_e[,j] * rowSums(z[k,]^2) + C[k] * tau_x[,j] + tau_w[1,j])))  ## D by K matrix
      for(k in 1:K){
        tau_wtk <- tau_wt[,k] ## vector of length D
        for(d in 1:D){
          res_d <- sapply(1:Gammai, function(j) sum(Y[d,] - (alpha[1,j] * w_tildet0[d,] + (1-alpha[1,j]) * wt_new[d,]) %*% z) ## vector of length Gammai
          mu_wt_dk <- sum(sapply(1:Gammai, function(j){
            tau_e[d,j] * (1 - alpha[1,j]) * sum(z[k,] * (res_d[j] + (1 - alpha[1,j]) * z[k,])) +
            tau_x[d,j] * sum(X[d,celltype == k]) + gamma[[j]][d,k] * v[[j]][d,k] * tau_w[1,j]))/tau_wtk[d]
          }
          wt_new[d,k] <- rnorm(mu_wt_dk, tau_wtk[d], 1)
        }
      }
    } else{
      wt_new <- w[[t]]
      tau_wt <- rowSums( Ck * tau_x + t(matrix(rep(tau_w, nrow(tau_x)), nrow=Gammai))  ## d by 1 matrix
      for(k in 1:K){
        mu_wtk <- rowSums(sapply(1:Gammai, function(j) tau_x[,j] * rowSums(X[,celltype == k]) + gamma[[j]][,k] * v[[j]][,k] * tau_w[1,j]))/tau_wt
        wt_new[,k] <- sapply(1:D, function(d) rnorm(mu_wtk[d], tau_wt[d], 1))
      }
    }
  }
  theta1_new <- list(w=w, z=z)
  return(list(theta1=theta1_new, theta2=theta2))
}
