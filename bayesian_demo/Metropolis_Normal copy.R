set.seed(532543)

#Simulating Data
theta_true<-14.65
sigma2<-1
mu<-0
tau2<-10000
n<-100
y<-rnorm(n,mean=theta_true,sd=sqrt(sigma2))

#Function to Calculate h(theta)
log_h<-function(theta_vec){
       vec<-rep(0,times=chains)
       for(j in 1:chains){
          vec[j]<- -(1/(2*sigma2))*sum((y-theta_vec[j])^2) - (1/(2*tau2))*((theta_vec[j]-mu)^2)
          }
       return(vec)
       }

mcmc_samples<-10000
chains<-3
theta<-matrix(0,nrow=mcmc_samples,ncol=chains)

#Starting Values
theta[1,]<-rnorm(n=chains,sd=100)

accept<-rep(1,times=chains)
#metrop_var<-10.00
metrop_var<-0.20
#metrop_var<-0.0001

for(i in 2:mcmc_samples){
   
   theta_proposed<-rnorm(n=chains,mean=theta[(i-1),],sd=sqrt(metrop_var))
   theta[i,]<-theta[(i-1),]
   ratio<-exp(log_h(theta_proposed) - log_h(theta[(i-1),]))
   for(j in 1:chains){
      if(ratio[j]>=runif(n=1,min=0,max=1)){
        theta[i,j]<-theta_proposed[j]
        accept[j]<-accept[j]+1
        }
      }
   
   print(accept/i)
   }

plot(theta[,1],type="l",ylim=c(-300,300))
for(j in 2:chains){
   lines(theta[,j],col=j)
   }


library(coda)

#Burn-in Amount
burnin<-2000
plot(beta[(burnin+1):mcmc_samples,1],type="l")
for(j in 2:chains){
   lines(beta[(burnin+1):mcmc_samples,j],col=j)
   }

#Autocorrelation Function Plot for Each Chain
par(mfrow=c(3,1))
for(j in 1:chains){
   acf(beta[(burnin+1):mcmc_samples,j])
   }

#Effective Sample Size for Each Chain
for(j in 1:chains){
   print(effectiveSize(beta[(burnin+1):mcmc_samples,j]))
   }

#Geweke Diagnostic
for(j in 1:chains){
   print(geweke.diag(beta[(burnin+1):mcmc_samples,j]))
   } 

#Gelman and Rubin Diagnostic
burnin<-2000
total_count<-(mcmc_samples-burnin)
keep_set1<-seq((burnin+1),burnin+(total_count/2),1)
keep_set2<-seq(burnin+(total_count/2)+1,mcmc_samples,1)

mh.list1<-
list(as.mcmc(theta[keep_set1,1]),
     as.mcmc(theta[keep_set2,1]),
     as.mcmc(theta[keep_set1,2]),
     as.mcmc(theta[keep_set2,2]),
     as.mcmc(theta[keep_set1,3]),
     as.mcmc(theta[keep_set2,3]))

mh.list <- mcmc.list(mh.list1)
gelman.diag(mh.list)$psrf
     
#Final Inference
full_beta<-
c(beta[(burnin+1):mcmc_samples,1],
  beta[(burnin+1):mcmc_samples,2],
  beta[(burnin+1):mcmc_samples,3])
  
print('full_beta')
print(head(full_beta))

saveRDS(full_beta, paste0('/home/wd262/scratch60/cytof/write/11272019/mcmc_beta_HL_',
          func_marker_to_be_tested,"_",celltype_test,'.rds'))

print('save full_beta')
                         
full_alpha<-list()
full_alpha[[1]]<-
c(alpha[[1]][(burnin+1):mcmc_samples,1],
  alpha[[1]][(burnin+1):mcmc_samples,2],
  alpha[[1]][(burnin+1):mcmc_samples,3])


full_alpha[[2]]<-
c(alpha[[2]][(burnin+1):mcmc_samples,1],
  alpha[[2]][(burnin+1):mcmc_samples,2],
  alpha[[2]][(burnin+1):mcmc_samples,3])



full_alpha[[3]]<-
c(alpha[[3]][(burnin+1):mcmc_samples,1],
  alpha[[3]][(burnin+1):mcmc_samples,2],
  alpha[[3]][(burnin+1):mcmc_samples,3])



full_alpha[[4]]<-
c(alpha[[4]][(burnin+1):mcmc_samples,1],
  alpha[[4]][(burnin+1):mcmc_samples,2],
  alpha[[4]][(burnin+1):mcmc_samples,3])

saveRDS(full_alpha, paste0('/home/wd262/scratch60/cytof/write/11272019/mcmc_alpha_HL_',
          func_marker_to_be_tested,"_",celltype_test,'.rds'))

                         
full_beta0<-
c(beta0[(burnin+1):mcmc_samples,1],
  beta0[(burnin+1):mcmc_samples,2],
  beta0[(burnin+1):mcmc_samples,3])

saveRDS(full_beta0, paste0('/home/wd262/scratch60/cytof/write/11272019/mcmc_beta0_HL_',
        func_marker_to_be_tested,"_",celltype_test,'.rds'))
                         
full_beta0_t<-
c(beta0_t[(burnin+1):mcmc_samples,1],
  beta0_t[(burnin+1):mcmc_samples,2],
  beta0_t[(burnin+1):mcmc_samples,3])

saveRDS(full_beta0_t, paste0('/home/wd262/scratch60/cytof/write/11272019/mcmc_beta0_t_HL_',
          func_marker_to_be_tested,"_",celltype_test,'.rds'))  

# Save W
saveRDS(data$W,paste0('/home/wd262/scratch60/cytof/write/11272019/mcmc_W_HL_',
        func_marker_to_be_tested,"_",celltype_test,'.rds'))