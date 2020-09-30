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

  tau_z <- lapply(1:Lambdai, function(j) colSums(tau_e[,j] * w_pseudo[[j]]^2) + 1 )
  tau_z <- Reduce('+', tau_z) # 1*K
  tau_z <- matrix(tau_z, nrow = K, ncol = 1)[,rep(1,N)]
  

  res <- lapply(w_pseudo, function(w_p) Y0 - w_p %*% z_new )
  mu_z <- lapply(1:Lambdai, function(j) t(tau_e[,rep(j,K)] * w_pseudo[[j]]) %*% res[[j]]  )
  mu_z <- Reduce('+', mu_z)/tau_z

  for(k in 1:K){
    for(n in 1:N){
      z_new[k,n] <- rnorm(1, mu_z[k,n], sd = sqrt(1/tau_z[k,n]))
    }
  }
  

  return(z_new)
  
}

update_w <- function(Y0, X0, W_tilde,
                      tau_e, alpha, gamma, v, tau_x, tau_w, z, w, t0=1
                      ){
  Lambdai = ncol(tau_e)
  end_idx = 0
  w_new <- list()
  for(t in 1:T){
    wt_new <- w[[t]]
    Xt <- X[[t]]
    ## here below w is updated by each element w_t0[d,k]
    C = as.vector(c_k[,t])
    start_idx <- end_idx + 1
    end_idx <- start_idx + sum(c_k[,t])-1

    tau_wt <- lapply(1:Lambdai, function(j) matrix(tau_x[,j], nrow = D, ncol = 1) %*% matrix(C, nrow =1, ncol=K)+tau_w[1,j])
    tau_wt <- Reduce('+', tau_wt)

    res1 <- sapply(c(1:K), function(k) apply(Xt[, Cl[start_idx:end_idx] == k],1,sum) )
    res1 <- tau_x[,rep(1,K)] * res1
    res2 <- lapply(1:Lambdai, function(j) tau_w[1,j]*v[[j]]*gamma[[j]])

    
    if(t == t0){
        tau_wt0 <- lapply(c(1:Lambdai), function(j) 
        (1-alpha[1,j])^2 * tau_e[,j] %*% (matrix(apply(z, 1,sum), nrow = 1, ncol = K))^2 )  ## Lambdai*D*K
        tau_wt <- tau_wt + Reduce('+', tau_wt0)
       
        #tau_wtk <- tau_wt[,k] ## vector of length D
        res <- lapply(1:Lambdai, function(j) Y0- alpha[1,j] * W_tilde  %*% z - (alpha[1,j] * W_tilde + (1-alpha[1,j]) * wt_new) %*% z)  ## list of number Lamdai of matrix of D*N, Lambdai*D*N
        res <- lapply(res, function(r) r %*% t(z)) # Lambdai*D*K

        mu_wt <- Reduce('+', res) + Reduce('+', res1) + Reduce('+', res2) # D*K
    } else{      
      mu_wt <- Reduce('+', res1) + Reduce('+', res2) # D*K
    }

    
    
    mu_wt <- mu_wt/tau_wt
    for(k in 1:K){
        for(d in 1:D){
            wt_new[d,k] <- rnorm(1, mu_wt[d,k], sd = sqrt(1/tau_wt[d,k]))
        }
    }

    w_new[[t]] <- wt_new

  }
  return(w_new)
}


