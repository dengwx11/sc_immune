mydir <- "/Users/wenxuandeng/GoogleDrive/sucksalt/SC/sc_immune/write/corrupt/"
rst <- readRDS(paste0(mydir,"rst.alpha=1.T=5.D=500.K=8.corrupt=0.01.tauW=1.tauXdBeta=1.rds"))

gene.name <- paste0("G",seq(500))
celltype.name <- paste0("C",seq(8))
sample.name <- paste0("S",seq(200))

mixture.demo <- rst$same_input$Y0
sigmatrix.demo <- rst$same_input$W_tilde

rownames(mixture.demo) <-  gene.name
rownames(sigmatrix.demo) <- gene.name
colnames(mixture.demo) <- sample.name
colnames(sigmatrix.demo) <- celltype.name


write.csv(mixture.demo, file =paste0(mydir, "mixture.demo.txt"),quote=FALSE)
write.csv(sigmatrix.demo, file =paste0(mydir, "sigmatrix.demo.txt"), quote=FALSE)




Docker run -v  absolute/path/to/input/dir:/Users/wenxuandeng/GoogleDrive/sucksalt/SC/sc_immune/write/corrupt -v absolute/path/to/output/dir:/Users/wenxuandeng/GoogleDrive/sucksalt/SC/sc_immune/write/cibersortx_output/corrupt --username wenxuan.deng@yale.edu --token b6dcbe0ff7466a84cd37f45662950113 --sigmatrix sigmatrix.demo.txt --mixture mixture.demo.txt --fraction 0 --rmbatchSmode FALSE