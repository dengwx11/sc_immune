# Simulation functions
#
# Simulate signature matrix W=D*K
# Simulate single cell read counts X=
# pi_ber=0.3
# alpha = 0.5
###################################################
##generate pi
#generate_pi <- function(alpha=0.001, beta=0.001){
#    pi_ber <- rbeta(n=1, alpha, beta)
#    return(pi_ber)
#}

#generate v
generate_v <- function(D, K, mu=2, tau=4){
  v <- matrix(abs(rnorm(D*K, mu, sd = (1/sqrt(tau)))), nrow = D, ncol = K)
  return(v)
}

#generate gamma
generate_gamma <- function(pi_ber, D, K
){
  gamma <- matrix(rbinom(n=D*K, size = 1, prob = pi_ber), nrow = D, ncol = K)
  return(gamma)
}

#generate w
generate_w <- function(v, gamma,D, K, tau_w=tau_w_para
){   
  w <- matrix(0, nrow = D, ncol = K)
  for (d in 1:D){
    for (k in 1:K){
      w[d,k] = max(rnorm(1, v[d,k]*gamma[d,k], sd = (1/sqrt(tau_w))),0.001)
    }
  }
  #w_new <- w*gamma
  #    print(paste("empirical pi=", (1-sum(w_new==0)/length(w_new))))
  return(w)
}

#calculate the max cos among cell types
max_cos <- function(w){
  K <- ncol(w)
  w_cos <- matrix(0, ncol= K, nrow = K)
  for (i in 1:K){
    for (j in 1:K){
      if (i==j){
        w_cos[i,j] = -1
      } else {
        w_cos[i,j] = (sum(w[,i]*w[,j]))/(sqrt(sum((w[,i]^2)))*sqrt(sum((w[,j]^2))))
      }
    }
  }
  max_cos = max(w_cos)
  return(max_cos)
}

#simulate w
simulate_sigmat_w <- function(iteration, D, K, pi_ber, T){
  w.list <- list()
  v.list <- list()
  gamma.list <- list()
  cos.list <- list()
  output <- list()
  for (i in 1:iteration){
    v <- generate_v(D, K)
    gamma <- generate_gamma(pi_ber, D, K)
    w <- generate_w(v, gamma, D, K)
    v.list[[i]] <- v
    gamma.list[[i]] <- gamma
    w.list[[i]] <- w
    cos.list[[i]] <- max_cos(w)
  }
  minindex <- which.min(cos.list)
  print(cos.list[[minindex]])
  # output$w <- w.list[[minindex]]
  output$v <- v.list[[minindex]]
  output$gamma <- gamma.list[[minindex]]
  output$max_cos <- cos.list[[minindex]]
  
  w.TissueList <- list()
  w.TissueList[[1]] <- w.list[[minindex]]
  if(T>1){
    for(i in 2:T) w.TissueList[[i]] <- generate_w(output$v, output$gamma, D, K)
  }
  output$w <- w.TissueList
  return(output)
}

#generate z
generate_z <- function(K,N){
  Z <- matrix(0, nrow = K, ncol = N)
  for (n in 1:N){
    zn = 0
    for (k in 1:(K-1)){
      Z[k,n] = runif(1, min=0, max=(1-zn))
      zn = zn + Z[k,n]
    }
    Z[K,n] = 1 - zn
  }
  return(Z)
}


#if wdkt>0, xdkc > 0
#if wdkt=0, mu=0, tau_xd=5, if x <0 , x =0.001
#generate X, a list contain read count matrix, cell labels, tau_xd, c_k, w_tilde, and celltype list
simulate_X <- function(D, K, w_output, T,corrupt_pi){
  
  X <- list()
  tau_xd = rgamma(D, 1, tau_xd_beta_para)
  gamma <- w_output$gamma
  for(t in 1:T){
    
    X[[t]] <- list()
    celltype.list <- c()
    w <- as.matrix(w_output$w[[t]])
    c_k <- sample(500:1000, size = K) #set the size for each cell type
    X[[t]]$counts <- matrix(0, nrow = D, ncol = sum(c_k))
    X[[t]]$Celltype_used <- rep("unknown", sum(c_k))
    X[[t]]$tau_xd <- tau_xd #sample random tau_xd from gamma distribution
    X[[t]]$c_k <- c_k
    X[[t]]$w_tilde <- matrix(0, nrow = D, ncol = K)
    end_idx = 0
    for (k in 1:K){
      celltype.list[k] <- paste("celltype",k)
      start_idx = end_idx + 1
      end_idx = start_idx + c_k[k] - 1
      for (d in 1:D){
        # if(t==1){
        #     if(gamma[d,k]==1){
        #       w[d,k]=w[d,k]*rbinom(1,1,1-corrupt_pi)
        #     }             
        #   #print(c(w_prev,w[d,k]))
        # }
        corrupt_flg <- rbinom(1,1,1-corrupt_pi)
        for (i in start_idx:end_idx){ ### please check here, if wdkt>0, xdkc < 0, should x_dkc = 0.001? or abs(rnorm()), or 0.
          X[[t]]$counts[d,i] = max(rnorm(1, w[d,k], sd = (1/sqrt(X[[t]]$tau_xd[d]))), 0.001 )
          if(gamma[d,k]==1){
            X[[t]]$counts[d,i] = X[[t]]$counts[d,i] * corrupt_flg
          }
        }
        X[[t]]$w_tilde[d,k] = mean(X[[t]]$counts[d,start_idx:end_idx])  ## observation
      }
      X[[t]]$Celltype_used[start_idx:end_idx] <- celltype.list[k]
    }
    X[[t]]$Celltype_list <- celltype.list
  }
  
  
  return(X)
}

corrupt_X <- function(X_sim_output, corrupt_pi = 0.2){    
  for(t in 1:T){
    X <- X_sim_output[[t]]
    #corrupt_bin <- matrix(rbinom(sum(X$c_k), 1, 1-corrupt_pi), nrow = 1, ncol = sum(X$c_k))[rep(1,D),]
    corrupt_bin <- matrix(rbinom(sum(X$c_k)*D, 1, 1-corrupt_pi), nrow = D, ncol = sum(X$c_k))
    X$counts <- X$counts * corrupt_bin
    c_k <- X$c_k
    end_idx = 0
    for (k in 1:K){
      start_idx = end_idx + 1
      end_idx = start_idx + c_k[k] - 1
      for (d in 1:D){                
        X$w_tilde[d,k] = mean(X$counts[d,start_idx:end_idx])
      }
    }
    
    X_sim_output[[t]] <- X
  }
  return(X_sim_output)
}


#X_sim_output <- corrupt_X(X_sim_output_original, corrupt_pi = 0.2)
#cbind(X_sim_output[[1]]$w_tilde[,2],X_sim_output_original[[1]]$w_tilde[,2])

#calculate W_tilde given X, if X does not contain w_tilde
#calculate_w_tilde <- function(X, D, K){
#    celltype.list <- X$Celltype_list
#    w_tilde <- matrix(0, nrow = D, ncol = K)
#    for (k in 1:k){
#        idx <- which(X$Celltype_used == celltype.list[k])
#        counts <- X$counts[,idx]
#        for (d in 1:D){
#            w_tilde[d,k] = mean(counts[d,])
#        }
#    }
#    return(w_tilde)
#}
#Y=W*Z+e 0.001 e mu_e=0, tau_e 
#simulate Y
simulate_y <- function(w_sim_output, X_sim_output, Z, D, N,
                       alpha = 0.5, mu_e = 0, tau_e = 1000){
  #W = (1 - alpha)*(w_sim_output$w[[1]]) + alpha*(X_sim_output$w_tilde)
  #W = w_sim_output$w[[1]]
  W = w_sim_output$v * w_sim_output$gamma
  e <- matrix(rnorm(D*N, mu_e, sd = (1/sqrt(tau_e))), nrow = D, ncol = N)
  e[e < 0] = 0.001
  Y = W %*% Z + e
  return(Y)
}
