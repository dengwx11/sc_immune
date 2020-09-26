source("SAME/Simulation.R")
set.seed(2020)

T=4
D=500
K=3
pi_ber = 0.3
N = 200 # bulk Y sample size
w_sim_output <- simulate_sigmat_w(500, D, K, pi_ber, T)
Z <- generate_z(K, N)
T=4
X_sim_output <- simulate_X(D, K, w_sim_output, T)
X<-list()
for(i in 1:T) X[[i]] <- X_sim_output[[i]]$counts
Y <- simulate_y(w_sim_output, X_sim_output[[1]], Z, D, N)

same_input <- list()
same_input$Y0 = Y
same_input$X = X
same_input$T = T
same_input$K = K
same_input$D = D


c_k = Reduce(function(d1,d2){
if(!is.numeric(d1)) return(cbind(as.matrix(d1$c_k,ncol=1),as.matrix(d2$c_k,ncol=1)))
else return(cbind(d1,as.matrix(d2$c_k,ncol=1)))
} ,X_sim_output)
same_input$c_k = c_k


C0 = sum(c_k)
same_input$W_tilde = X_sim_output[[i]]$w_tilde
#same_input$YSG 