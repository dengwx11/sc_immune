set.seed(2021)

update_w <- function(X_mat, 
                      gamma, v, tau_x, tau_w, w, Cl
                      ){
  Lambdai = ncol(tau_w)
  end_idx = 0
  w_new <- list()
  
  for(t in 1:T){
    wt_new <- w[[t]]
    Xt <- X_mat[[t]]
    ## here below w is updated by each element w_t0[d,k]
    C = as.vector(c_k[,t])
    
    start_idx <- end_idx + 1
    end_idx <- start_idx + sum(c_k[,t])-1

    #C_mat <- matrix(C, nrow =1, ncol=K)
    sample.cnt <- 500
    C_mat <- matrix(sample.cnt, nrow =1, ncol=K)

    
    tau_wt <- lapply(1:Lambdai, function(j) matrix(tau_x[,j], nrow = D, ncol = 1) %*% C_mat +tau_w[1,j])
    tau_wt <- Reduce('+', tau_wt)

    k.idx <- which(Cl[start_idx:end_idx] == k)

    if(length(k.idx)>0){
      res1 <- sapply(c(1:K), function(k) apply(Xt[, sample(k.idx,sample.cnt,replace=T)],1,sum))

      res1 <- lapply(1:Lambdai, function(j) tau_x[,rep(j,K)] * res1)
    }else{
      res1 <- lapply(1:Lambdai, function(j) matrix(0,nrow=D,ncol=K))
    }


    res2 <- lapply(1:Lambdai, function(j) tau_w[1,j]*v[[j]]*gamma[[j]])

    
      
    mu_wt <- Reduce('+', res1) + Reduce('+', res2) # D*K


    
    
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
