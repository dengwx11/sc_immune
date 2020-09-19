source("Simulation.R")
set.seed(2020)

w_sim_output <- simulate_sigmat_w(500, 500, 3, 0.3)
Z <- generate_z(3, 200)
X_sim_output <- simulate_X(500, 3, w_sim_output)
Y <- simulate_y(w_sim_output, X_sim_output, z, 500, 200)