corr_celltype <- c()
pr_gamma_rst <- c()
empirical_pi <- c(0.1, 0.15, 0.2, 0.25, 0.3)
rst_new.list <- list()
for (i in 1:length(empirical_pi)) {
  rst_counts_alpha1 <- SAME(Y0_sc, X, W_tilde, empirical_pi[i], mcmc_samples_theta1, Lambda, c_k, YSG, alpha =1) #time  554.109s
  
  #calculate v*g
  est_vg <- lapply(c(1:mcmc_samples_theta1), function(i) rst_counts_alpha1$theta2$gamma[[i]] * rst_counts_alpha1$theta2$v[[i]])
  average_est_vg_1 <- Reduce('+',est_vg)/mcmc_samples_theta1
  rst_counts_alpha1$theta2$vg <- average_est_vg_1
  #estimate z by nnls
  z_est_vg_nnls <- t(sapply(c(1:N), function(j) nnls(average_est_vg_1,(Y0_sc[YSG,j]))$x))
  corr_vg_nnls<- sapply(c(1:4), function(i) cor(as.vector(as.matrix(gt_sub)[,i]), as.vector(z_est_vg_nnls[,i])))
  corr_vg_all <- cor(as.vector(as.matrix(gt_sub)), as.vector(z_est_vg_nnls))
  corr_all <- append(corr_all, corr_vg_all)
  #corr_celltype[i,] <- corr_vg_nnls
  rst_counts_alpha1$theta2$corr_celltype <- corr_vg_nnls
  rst_counts_alpha1$theta2$corr_all <- corr_vg_all
  rst_new.list[[i]] <- rst_counts_alpha1
  saveRDS(rst_counts_alpha1, file = paste0("SAME/V2/rst_tauv0.1_updatepi", empirical_pi[i], "_alpha1.rds"))
}
names(rst_new.list) <- paste0("pi",empirical_pi)

gamma_tmp <- lapply(1:5, function(j) unlist(lapply(rst_new.list[[j]]$theta2$gamma, function(x) sum(x))))
for (i in 1:5){
  corr_celltype <- c(corr_celltype,  rst_new.list[[i]]$theta2$corr_celltype)
  pr_gamma <- gamma_tmp[[i]]/(D*K)
  pr_gamma_rst <- c(pr_gamma_rst, pr_gamma)
}
corr_celltype <- c(corr_celltype, corr_wtilde)

df_corr_ct <- data.frame(corr = corr_celltype,
                         celltype = rep(celltype.list, 60),
                         initial_pi = rep(c(rep(rep(c(names(rst_new.list), "W_tilde"), each = 4), 5)), 2),
                         tau_v = rep(c(rep(0, 24), rep(1, 24), rep(0.1, 24), rep(4, 24), rep(0.01, 24)),2),
                         update_pi = c(rep("original", 120), rep("noupdate",120))
)
df_corr_ct$update_pi_tau_v <- paste0(df_corr_ct$update_pi, "_", df_corr_ct$tau_v)
saveRDS(df_corr_ct, file = "SAME/V2/corr_compare_df.rds")


p <- ggplot(df_corr_ct, aes(x=initial_pi, y=corr, colour = celltype)) + geom_jitter(width = 0.1) + facet_wrap(~update_pi_tau_v, nrow =2)

p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="red")
p1 <- ggplot(data=df_corr_ct,aes(x=initial_pi, y=corr, group=celltype, colour = celltype))+geom_line()

iteration = rep(c("1-10", "11-20", "21-30", "31-40", "41-50"), each = 10)
df_gamma <- data.frame(Pr = pr_gamma_rst,
                       initial_pi = rep(rep(rep(names(rst_new.list), each = 50), 5),2),
                       iteration = rep(rep(rep(iteration, 5),5),2),
                       tau_v = rep(c(rep(0, 250), rep(1, 250), rep(0.1, 250), rep(4, 250), rep(0.01, 250)),2),
                       update_pi = c(rep("original", 1250), rep("noupdate",1250))
)
df_gamma$update_pi_tau_v <- paste0(df_gamma$update_pi, "_", df_gamma$tau_v)
p2 <- ggplot(df_gamma, aes(x=initial_pi, y=Pr, colour = iteration)) + geom_jitter(width = 0.1)+ facet_wrap(~update_pi_tau_v, nrow = 2)
p2 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red")
saveRDS(df_gamma, file = "SAME/V2/gamma_compare_df.rds")
df_gamma <- readRDS(file = "SAME/V2/gamma_compare_df.rds")