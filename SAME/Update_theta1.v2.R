set.seed(2021)

update_w <- function(Y0, X_mat, W_tilde,
                      tau_e, alpha, gamma, v, tau_x, tau_w, z, w, t0=1, Cl
                      ){
  Lambdai = ncol(tau_e)
  end_idx = 0
  w_new <- list()
  
  for(t in 1:T){
    wt_new <- w[[t]]
    Xt <- X_mat[[t]]
    ## here below w is updated by each element w_t0[d,k]
    C = as.vector(c_k[,t])
    
    start_idx <- end_idx + 1
    end_idx <- start_idx + sum(c_k[,t])-1
    
    tau_wt <- lapply(1:Lambdai, function(j) matrix(tau_x[,j], nrow = D, ncol = 1) %*% matrix(C, nrow =1, ncol=K)+tau_w[1,j])
    tau_wt <- Reduce('+', tau_wt)

    res1 <- sapply(c(1:K), function(k) apply(Xt[, Cl[start_idx:end_idx] == k],1,sum) )

    res1 <- lapply(1:Lambdai, function(j) tau_x[,rep(j,K)] * res1)

    res2 <- lapply(1:Lambdai, function(j) tau_w[1,j]*v[[j]]*gamma[[j]])

    
    if(t == t0){
        tau_wt0 <- lapply(c(1:Lambdai), function(j) 
        (1-alpha[1,j])^2 * tau_e[,j] %*% matrix(apply(z^2, 1,sum), nrow = 1, ncol = K) )  ## Lambdai*D*K
        tau_wt <- tau_wt + Reduce('+', tau_wt0)
       
        #tau_wtk <- tau_wt[,k] ## vector of length D
        res <- lapply(1:Lambdai, function(j) Y0- (alpha[1,j] * W_tilde + (1-alpha[1,j]) * wt_new) %*% z)  ## list of number Lamdai of matrix of D*N, Lambdai*D*N
        res <- lapply(res, function(r) r %*% t(z)) # Lambdai*D*K
        res0 <- lapply(1:Lambdai, function(j) (1-alpha[1,j]) * wt_new * matrix(apply(z^2,1,sum), nrow = 1, ncol=K)[rep(1,D),]) # Lambdai*D*K
        res <- lapply(1:Lambdai, function(j) res[[j]]+res0[[j]])
        res <- lapply(1:Lambdai, function(j) res[[j]] * (1-alpha[1,j]) * tau_e[,rep(j,K)])


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

    wt_new <- wt_new * (wt_new>0)

    w_new[[t]] <- wt_new

  }
  return(w_new)
}
