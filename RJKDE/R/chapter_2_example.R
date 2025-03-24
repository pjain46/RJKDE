rnormmix <- function(n, mu, sigma,w){
  k <- length(mu)
  inds <- sample(k, n, prob = w, replace=TRUE)
  
  xsamp <- rnorm(n,mu[inds], sigma[inds])
  return(list(inds = inds, xsamp = xsamp))
}

normmix_density <- function(x, mus, sigs, q){
  
  A <- sum(q*dnorm(x, mus, sigs))
  return(A)
}


library(dirmult)
library(rnn)
library(Rcpp)

Rcpp::sourceCpp("~/Dropbox (ASU)/Palak_Hahn/RJMIXMOD.cpp")


library(rnn)

# 2021 Data

year <- 2021
file_path <- sprintf("~/Dropbox (ASU)/Palak_Hahn/Birthweight Data/natality%dus.csv", year)
natalitydata <- read.csv(file_path)
names(natalitydata)


# BOY: sex
# MARRIED: dmar
# BLACK: mrace15 (=2 means Black Only)
# OVER33: mager
# HIGH SCHOOL: meduc
# FULL PRENATAL: precare5
# SMOKER: cig_0
# BIRTH WEIGHT: dbwt

dat <- natalitydata[,c('sex', 'dmar', 'mrace15', 'mager', 'meduc', 'precare5', 'cig_0', 'dbwt')]

dat <- na.omit(dat)
y <- log(dat$dbwt/1000)

# xnode <- apply(X[, 1:p], 1, function(a) bin2int(t(as.matrix(as.numeric(a), 1, p))) + 1)

# sdy <- sd(y)
# y <- y/sdy

gridsize <- 500
ygrid <- seq(min(y),max(y),length.out = gridsize)




mc <- 5000
kinit <- 100
burn <- 1
kde <- density(sample(y,kinit))
bw <- kde$bw

print(bw)

KDEpost <- rj_mcmc_rcpp(y, ygrid, k = kinit, bw = bw, mc = mc, mu_step_size = 0.1)
fsamps <- KDEpost$fsamps
fsamps <- fsamps[,(burn+1):mc]

fhat <- rowMeans(fsamps)
fhat.u <- apply(fsamps,1,quantile,0.975)
fhat.l <- apply(fsamps,1,quantile,0.025)

ygrid <- KDEpost$ygrid
hist(y,freq=FALSE,35,main='',xlab='y')
lines(ygrid,fhat,col='green',lwd=2,lty=2)
lines(ygrid,fhat.u,lty=3,col='red')
lines(ygrid,fhat.l,lty=3,col='red')


fsamps <- KDEpost$fsamps/exp(ygrid)
fsamps <- fsamps[,(burn+1):mc]

fhat <- rowMeans(fsamps)
fhat.u <- apply(fsamps,1,quantile,0.975)
fhat.l <- apply(fsamps,1,quantile,0.025)

ygrid <- KDEpost$ygrid
hist(exp(y),freq=FALSE,35,main='',xlab='Y (in Kgs)', xlim=c(0,5))
lines(exp(ygrid),fhat,col='green',lwd=2,lty=2)
lines(exp(ygrid),fhat.u,lty=3,col='red')
lines(exp(ygrid),fhat.l,lty=3,col='red')






# 
# 
# fsamps <- exp(KDEpost$fsamps)
# fsamps <- (fsamps - min(fsamps))/(max(fsamps) - min(fsamps))
# 
# fhat <- rowMeans(fsamps)
# fhat.u <- apply(fsamps,1,quantile,0.975)
# fhat.l <- apply(fsamps,1,quantile,0.025)
# 
# 
# hist(exp(y*sdy),ylim = c(0,1),freq=FALSE,35,main='',xlab='y')
# lines(exp(ygrid*sdy),fhat,col='green',lwd=2,lty=2)
# lines(exp(ygrid*sdy),fhat.u,lty=3,col='red')
# lines(exp(ygrid*sdy),fhat.l,lty=3,col='red')
# 
# 
# 



