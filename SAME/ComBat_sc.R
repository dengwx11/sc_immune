#' Adjust for batch effects using an empirical Bayes framework
#' 
#' ComBat allows users to adjust for batch effects in datasets where the batch covariate is known, using methodology 
#' described in Johnson et al. 2007. It uses either parametric or non-parametric empirical Bayes frameworks for adjusting data for 
#' batch effects.  Users are returned an expression matrix that has been corrected for batch effects. The input
#' data are assumed to be cleaned and normalized before batch effect removal. 
#' 
#' @param dat Genomic measure matrix (dimensions probe x sample) - for example, expression matrix
#' @param batch {Batch covariate (only one batch allowed)}
#' @param to.sc TRUE If FALSE ComBat_sc is the same function with ComBat
#' @param mod Model matrix for outcome of interest and other covariates besides batch
#' @param par.prior (Optional) TRUE indicates parametric adjustments will be used, FALSE indicates non-parametric adjustments will be used
#' @param prior.plots (Optional) TRUE give prior plots with black as a kernel estimate of the empirical batch effect density and red as the parametric
#' @param mean.only (Optional) FALSE If TRUE ComBat only corrects the mean of the batch effect (no scale adjustment)
#'
#' @return If to.sc is TRUE, ComBat_sc returns a list composed of a bulk seq matrix adjusted to single cell space and empirical tau_e.
#' 
#' @export
#' 

ComBat_sc <- function(dat, batch, to.sc = TRUE, mod=NULL, par.prior=TRUE,prior.plots=FALSE,mean.only=FALSE) {
  # make batch a factor and make a set of indicators for batch
  dat <- as.matrix(dat)
  batch <- as.factor(batch)
  zero.rows.lst <- lapply(levels(batch), function(batch_level) {
    if (sum(batch == batch_level) > 1) {
      return(which(apply(dat[, batch == batch_level], 1, 
                         function(x) {
                           var(x) == 0
                         })))
    }
    else {
      return(which(rep(1, 3) == 2))
    }
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(dat), zero.rows)
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n", 
                length(zero.rows)))
    dat.orig <- dat
    dat <- dat[keep.rows, ]
  }
  if(mean.only==TRUE){cat("Using the 'mean only' version of ComBat\n")}
  if(length(dim(batch))>1){stop("This version of ComBat only allows one batch variable")}  ## to be updated soon!
  batch <- as.factor(batch)
  batchmod <- model.matrix(~-1+batch)  
  cat("Found",nlevels(batch),'batches\n')
  
  # A few other characteristics on the batches
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch){batches[[i]] <- which(batch == levels(batch)[i])} # list of samples in each batch  
  n.batches <- sapply(batches, length)
  if(any(n.batches==1)){mean.only=TRUE; cat("Note: one batch has only one sample, setting mean.only=TRUE\n")}
  n.array <- sum(n.batches)
  
  #combine batch variable and covariates
  design <- cbind(batchmod,mod)
  
  # check for intercept in covariates, and drop if present
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[,!check])
  
  # Number of covariates or covariate levels
  cat("Adjusting for",ncol(design)-ncol(batchmod),'covariate(s) or covariate level(s)\n')
  
  # Check if the design is confounded
  if(qr(design)$rank<ncol(design)){
    #if(ncol(design)<=(n.batch)){stop("Batch variables are redundant! Remove one or more of the batch variables so they are no longer confounded")}
    if(ncol(design)==(n.batch+1)){stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")}
    if(ncol(design)>(n.batch+1)){
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){stop('The covariates are confounded! Please remove one or more of the covariates so the design is not confounded')
      }else{stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")}}
  }
  
  ## Check for missing values
  NAs = any(is.na(dat))
  if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
  #print(dat[1:2,])
  ##Standardize Data across genes
  cat('Standardizing Data across genes\n')
  if (!NAs){B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))}else{B.hat=apply(dat,1,Beta.NA,design)} #Standarization Model
  grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
  if (!NAs){var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)}else{var.pooled <- apply(dat-t(design%*%B.hat),1,var,na.rm=T)}
  
  stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
  if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
  s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array))) #Z_ijg
    
  ##Get regression batch effect parameters
  cat("Fitting L/S model and finding priors\n")
  batch.design <- design[,1:n.batch]
  if (!NAs){
    gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))
  } else{
    gamma.hat=apply(s.data,1,Beta.NA,batch.design)
    
  }
  delta.hat <- NULL
  for (i in batches){
    if(mean.only==TRUE){delta.hat <- rbind(delta.hat,rep(1,nrow(s.data)))}else{
      delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
    }
  }
  
  ##Find Priors
  gamma.bar <- apply(gamma.hat, 1, mean)
  t2 <- apply(gamma.hat, 1, var)
  a.prior <- apply(delta.hat, 1, aprior)
  b.prior <- apply(delta.hat, 1, bprior)
  
  
  ##Plot empirical and parametric priors
  
  # if (prior.plots & par.prior){
  #   par(mfrow=c(2,2))
  #   tmp <- density(gamma.hat[1,])
  #   plot(tmp,  type='l', main="Density Plot")
  #   xx <- seq(min(tmp$x), max(tmp$x), length=100)
  #   lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
  #   qqnorm(gamma.hat[1,])	
  #   qqline(gamma.hat[1,], col=2)	
  #   
  #   tmp <- density(delta.hat[1,])
  #   invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
  #   tmp1 <- density(invgam)
  #   plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
  #   lines(tmp1, col=2)
  #   qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
  #   lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
  #   title('Q-Q Plot')
  # }
  
  ##Find EB batch adjustments
  
  gamma.star <- delta.star <- NULL
  if(par.prior){
    cat("Finding parametric adjustments\n")
    for (i in 1:n.batch){
      if(mean.only){
        gamma.star <- rbind(gamma.star,postmean(gamma.hat[i,],gamma.bar[i],1,1,t2[i]))
        delta.star <- rbind(delta.star,rep(1,nrow(s.data)))
      }else{temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
      }
    }
  }else{
    cat("Finding nonparametric adjustments\n")
    for (i in 1:n.batch){
      if(mean.only){delta.hat[i,]=1}
      temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
    }
  }
  
  
  ### Normalize the Data ###
  cat("Adjusting the Data\n")
  
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
    j <- j+1
  }
  
  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
  rst <- bayesdata
  
  #adjust the bulk seq data to single cell space
  if(to.sc == TRUE){
    cat("Adjusting the Data to single cell space\n")
    rst <- list()
    i <- batches[[1]]
    bayesdata <- ((s.data[,i]-t(batch.design[i,]%*%gamma.star))*(sqrt(delta.star[2,])%*%t(rep(1,n.batches[1]))))/(sqrt(delta.star[1,])%*%t(rep(1,n.batches[1])))
    stand.mean <- stand.mean[,i]
    i <- batches[[2]]
    bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.batches[1])))) + t(batch.design[i,]%*%gamma.star)[,1:n.batches[1]] + stand.mean
    rst$bulk_to_sc <- bayesdata
    rst$tau_e <- 1/(delta.star[2,]*var.pooled)
  }
  if (length(zero.rows) > 0) {
    dat.orig[keep.rows, 1:n.batches[1]] <- bayesdata
    bayesdata <- dat.orig[, 1:n.batches[1]]
    rst$bulk_to_sc <- bayesdata
  }
  
  return(rst)
  
}


sva.class2Model <- function(classes) {
  return(model.matrix(~factor(classes)))
}


modefunc <- function(x) {
  return(as.numeric(names(sort(-table(x)))[1]))
}

mono <- function(lfdr){.Call("monotone",as.numeric(lfdr),PACKAGE="sva")}


edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, lambda=0.8, ...) {
  pi0 <- mean(p >= lambda)/(1 - lambda)
  pi0 <- min(pi0, 1)
  
  n = length(p)
  transf = match.arg(transf)
  
  if(transf=="probit") {
    p = pmax(p, eps)
    p = pmin(p, 1-eps)
    x = qnorm(p)
    myd = density(x, adjust=adj)
    mys = smooth.spline(x=myd$x, y=myd$y)
    y = predict(mys, x)$y
    lfdr = pi0*dnorm(x)/y
  }
  
  if(transf=="logit") {
    x = log((p+eps)/(1-p+eps))
    myd = density(x, adjust=adj)
    mys = smooth.spline(x=myd$x, y=myd$y)
    y = predict(mys, x)$y
    dx = exp(x)/(1+exp(x))^2
    lfdr = pi0*dx/y
  }
  
  if(trunc) {lfdr[lfdr > 1] = 1}
  if(monotone) {	
    lfdr = lfdr[order(p)]
    lfdr = mono(lfdr)
    lfdr = lfdr[rank(p)]
  }
  
  return(lfdr)
}



# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
trim.dat <- function(dat){
  tmp <- strsplit(colnames(dat),'\\.')
  tr <- NULL
  for (i in 1:length(tmp)){tr <- c(tr,tmp[[i]][1]!='X')}
  tr
}

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2}
bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}


# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
  n <- apply(!is.na(sdat),1,sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while(change>conv){
    g.new <- postmean(g.hat,g.bar,n,d.old,t2)
    sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
    d.new <- postvar(sum2,n,a,b)
    change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count+1
  }
  #cat("This batch took", count, "iterations until convergence\n")
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star","d.star")
  adjust
}

#likelihood function used below
L <- function(x,g.hat,d.hat){prod(dnorm(x,g.hat,sqrt(d.hat)))}

# Monte Carlo integration functients
int.eprior <- function(sdat,g.hat,d.hat){
  g.star <- d.star <- NULL
  r <- nrow(sdat)
  for(i in 1:r){
    g <- g.hat[-i]
    d <- d.hat[-i]		
    x <- sdat[i,!is.na(sdat[i,])]
    n <- length(x)
    j <- numeric(n)+1
    dat <- matrix(as.numeric(x),length(g),n,byrow=T)
    resid2 <- (dat-g)^2
    sum2 <- resid2%*%j
    LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
    LH[LH=="NaN"]=0
    g.star <- c(g.star,sum(g*LH)/sum(LH))
    d.star <- c(d.star,sum(d*LH)/sum(LH))
    #if(i%%1000==0){cat(i,'\n')}
  }
  adjust <- rbind(g.star,d.star)
  rownames(adjust) <- c("g.star","d.star")
  adjust	
} 

#fits the L/S model in the presence of missing data values

Beta.NA = function(y,X){
  des=X[!is.na(y),]
  y1=y[!is.na(y)]
  B <- solve(t(des)%*%des)%*%t(des)%*%y1
  B
}

