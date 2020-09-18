# Simulation functions
#
# Simulate signature matrix W=D*K
# alpha=0.001
# beta=0.001
# pi_ber=0.5
#
###################################################
##generate pi
#generate_pi <- function(alpha=0.001, beta=0.001){
#    pi_ber <- rbeta(n=1, alpha, beta)
#    return(pi_ber)
#}

#generate v
generate_v <- function(D, K, mu=0, tau=0.01){
    v <- matrix(rnorm(D*K, mu, sd = (1/sqrt(tau))), nrow = D, ncol = K)
    return(v)
}

#generate gamma
generate_gamma <- function(pi_ber, D, K
                            ){
    gamma <- matrix(rbinom(n=D*K, size = 1, prob = pi_ber), nrow = D, ncol = K)
    return(gamma)
}

#generate w
generate_w <- function(v, gamma,D, K, tau_w=0.01
                            ){   
    w <- matrix(0, nrow = D, ncol = K)
#    tau_w <- rgamma(1, alpha, beta)
    for (d in 1:D){
        for (k in 1:K){
            w[d,k] = rnorm(1, v[d,k], sd = (1/sqrt(tau_w)))
        }
    }
    w_new <- w*gamma
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

simulation_function <- function(iteration, D, K, pi_ber){
    w.list <- list()
    cos.list <- list()
    for (i in 1:iteration){
        v <- generate_v(D, K)
        gamma <- generate_gamma(pi_ber, D, K)
        w <- generate_w(v, gamma, D, K)
        w.list[[i]] <- w
        cos.list[[i]] <- max_cos(w)
    }
    minindex <- which.min(cos.list)
    print(cos.list[minindex])
    w_simulated <- w.list[minindex]
    return(w_simulated)
}