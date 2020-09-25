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
Y <- simulate_y(w_sim_output, X_sim_output[[1]], Z, D, N)

same_input <- list()
same_input$Y0 = Y
same_input$X = X
same_input$T = T
same_input$K = K
same_input$D = D
c_k = same_input$c_k
C0 = sum(c_k)
W_tilde = same_input$W_tilde
YSG =  same_input$YSG