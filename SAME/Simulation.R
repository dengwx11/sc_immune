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
generate_v <- function(D, K, mu=1, tau=4){
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
generate_w <- function(v, gamma,D, K, tau_w=4
                            ){   
    w <- matrix(0, nrow = D, ncol = K)
    for (d in 1:D){
        for (k in 1:K){
            w[d,k] = rnorm(1, v[d,k], sd = (1/sqrt(tau_w)))
            if (w[d,k]<0){w[d,k]=0}
        }
    }
    w_new <- w*gamma
#    print(paste("empirical pi=", (1-sum(w_new==0)/length(w_new))))
    return(w_new)
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
simulate_sigmat_w <- function(iteration, D, K, pi_ber){
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
    output$w <- w.list[[minindex]]
    output$v <- v.list[[minindex]]
    output$gamma <- gamma.list[[minindex]]
    output$max_cos <- cos.list[[minindex]]
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
simulate_X <- function(D, K, w_output){
    celltype.list <- c()
    w <- as.matrix(w_output$w)
    c_k <- sample(0:200, size = K) #set the size for each cell type
    X <- list()
    X$counts <- matrix(0, nrow = D, ncol = sum(c_k))
    X$Celltype_used <- rep("unknown", sum(c_k))
    X$tau_xd <- rgamma(D, 1, 1) #sample random tau_xd from gamma distribution
    X$c_k <- c_k
    X$w_tilde <- matrix(0, nrow = D, ncol = K)
    end_idx = 0
    for (k in 1:K){
        celltype.list[k] <- paste("celltype",k)
        start_idx = end_idx + 1
        end_idx = start_idx + c_k[k] - 1
        for (d in 1:D){
            for (i in start_idx:end_idx){ ### please check here, if wdkt>0, xdkc < 0, should x_dkc = 0.001? or abs(rnorm()), or 0.
                if (w[d,k]==0){
                    X$counts[d,i] = rnorm(1, w[d,k], sd = (1/sqrt(X$tau_xd[d])))
                    if (X$counts[d,i]<0){X$counts[d,i] = 0.001}
                }else{
                    X$counts[d,i] = abs(rnorm(1, w[d,k], sd = (1/sqrt(X$tau_xd[d])))) ### there is rarely 0 in X
                }
            }
            X$w_tilde[d,k] = mean(X$counts[d,start_idx:end_idx])
        }
        X$Celltype_used[start_idx:end_idx] <- celltype.list[k]
    }
    X$Celltype_list <- celltype.list
    return(X)
}

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
                        alpha = 0.5, mu_e = 0, tau_e = 100){
    W = (1 - alpha)*(w_sim_output$w) + alpha*(X_sim_output$w_tilde)
    e <- matrix(rnorm(D*N, mu_e, sd = (1/sqrt(tau_e))), nrow = D, ncol = N)
    e[e < 0] = 0.001
    Y = W %*% Z + e
    return(Y)
}