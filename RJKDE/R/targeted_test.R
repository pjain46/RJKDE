library(Rcpp)
#sourceCpp("RJMIXMOD.cpp")



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

# generate data

n <- 50000

#mu = c(-1,0,2)
#sig = c(0.125,0.25,0.625)

kprime = 4
mu = 2*rbeta(kprime, 0.8,0.8)-0.5
sig = seq(0.125,0.75, length.out = kprime)


mu_h = -2.75
sig_h = 1.25

s = sqrt(1/(1/sig^2 - 1/sig_h^2))
m = (mu/sig^2 - mu_h/sig_h^2)*s^2

a = dnorm(mu_h,m,sqrt(s^2 + sig_h^2))

p = 1/a
p = p/sum(p)



ind <- sample(1:kprime,n,replace=TRUE,prob = p)

x = m[ind] + s[ind]*rnorm(n)

w <- dnorm(x,mu_h,sig_h,log=TRUE)
w <- w - max(w)
w <- exp(w)/sum(exp(w))

n.sub <- 300
y <- sample(x,n.sub,replace=TRUE,prob = w)
y2 <- sample(x,n.sub,replace=TRUE)





true_mu <- m
true_sig <- s
true_w <- p

rho <- rep(1,kprime)/kprime

sdy <- sd(y)
y <- y/sdy

sdy2 <- sd(y2)
y2 <- y2/sdy2






gridsize <- 500
burn <- 1000
ygrid <- seq(min(y)-1*sd(y),max(y)+1*sd(y),length.out = gridsize)
ygrid2 <- seq(min(y2)-1*sd(y2),max(y2)+1*sd(y2),length.out = gridsize)



ftemp <- sapply(ygrid*sdy,function(x) normmix_density(x, true_mu, true_sig,true_w))

ftemp_alt <- sapply(ygrid*sdy,function(x) normmix_density(x, mu, sig,rho))


ftemp2 <- sapply(ygrid2*sdy2,function(x) normmix_density(x, true_mu, true_sig,true_w))



hist(y*sdy,freq=FALSE,ylim=c(0,1), 35,main='',xlab='y')
lines(ygrid*sdy,ftemp,col='purple',lwd=2)
lines(ygrid*sdy,ftemp_alt,col='orange')
legend(1.5, 0.85, legend=c("True Density Curve "),col=c('purple'), lty=1,lwd = 2, cex=0.8)

plot(ygrid*sdy,ftemp,col='purple',lwd=2,xlab='y', ylab = 'Density', type = 'l')
legend(-1.5, 0.5, legend=c("True curve "),col=c('purple'), lty=1,lwd = 2, cex=0.8)

mc <- 20000
kinit = 20
kde = density(sample(y,kinit))
bw = kde$bw

print(bw)
print(sig_h/sdy)

kde2 = density(sample(y2,kinit))
bw2 = kde2$bw


  start_time <- Sys.time()
  
 # kinit = ceiling(sqrt(12*(1-(0.1*sig_h/sdy)^2) + 1))
  kinit = 10

  KDEpost <- rj_mcmc_rcpp(y, ygrid, k = kinit, bw = 0.5*sig_h/sdy, mc = mc, mu_prior_b = 10, mu_h = mu_h/sdy, sig_h = sig_h/sdy)

  KDEpost2 <- rj_mcmc_rcpp(y2, ygrid2, k = kinit, bw = 0.5, mc = mc, mu_prior_b = 10, mu_h = -100, sig_h = 1)
  
    
  print(Sys.time() - start_time)
  
  fsamps <- KDEpost$fsamps/sdy
  fsamps <- fsamps[,(burn+1):mc]
  
  
  fhat <- rowMeans(fsamps)
  fhat.u <- apply(fsamps,1,quantile,0.975)
  fhat.l <- apply(fsamps,1,quantile,0.025)
  
  fsamps_adj <- KDEpost$fsamps_adj/sdy
  fsamps_adj <- fsamps_adj[,(burn+1):mc]
  
  
  fhat_adj <- rowMeans(fsamps_adj,na.rm = TRUE)
  fhat_adj.u <- apply(fsamps_adj,1,quantile,0.975,na.rm = TRUE)
  fhat_adj.l <- apply(fsamps_adj,1,quantile,0.025,na.rm = TRUE)
  
  fsamps2 <- KDEpost2$fsamps/sdy2
  fsamps2 <- fsamps2[,(burn+1):mc]
  
  
  fhat2 <- rowMeans(fsamps2)
  fhat2.u <- apply(fsamps2,1,quantile,0.975)
  fhat2.l <- apply(fsamps2,1,quantile,0.025)
  
  
  
  
  par(mfrow=c(2,2), mai = c(1, 0.5, 0.05, 0.1))
  
  
  ygrid <- KDEpost$ygrid
  hist(y*sdy,freq=FALSE,ylim=c(0,1), 35,main='',ylab='y')
  lines(ygrid*sdy,ftemp,col='purple',lwd=2)
  lines(ygrid*sdy,ftemp_alt,col='orange')
  lines(ygrid2*sdy2,ftemp2,col='red',lwd=1,lty=2)
  lines(ygrid*sdy,fhat,col='green',lwd=2,lty=2)
  lines(ygrid*sdy,fhat.u,lty=3,col='red')
  lines(ygrid*sdy,fhat.l,lty=3,col='red')
  
  
    
  hist(y2*sdy2,freq=FALSE,ylim=c(0,1), 35,main='',ylab='y')
  ygrid2 <- KDEpost2$ygrid
  lines(ygrid2*sdy2,ftemp2,col='red',lwd=2)
  lines(ygrid2*sdy2,fhat2,col='green',lwd=2,lty=2)
  lines(ygrid2*sdy2,fhat2.u,lty=3,col='red')
  lines(ygrid2*sdy2,fhat2.l,lty=3,col='red')
  
  
  
  

  
  plot(ygrid[ygrid<0]*sdy,ftemp[ygrid<0],col='purple',lwd=4,type = 'l',xlim=c(-1,1)-1,bty='n',main='',ylab='y')
  lines(ygrid[ygrid<0]*sdy,fhat_adj[ygrid<0],col='green')
  lines(ygrid[ygrid<0]*sdy,fhat_adj.u[ygrid<0],lty=3,col='red',lwd=2)
  lines(ygrid[ygrid<0]*sdy,fhat_adj.l[ygrid<0],lty=3,col='red',lwd=2)
  
  
 
 
 
  
  plot(ygrid2[ygrid2<0]*sdy2,ftemp2[ygrid2<0],col='purple',lwd=4,type = 'l',xlim=c(-1,1)-1,bty='n',main='',ylab='y')
  lines(ygrid2[ygrid2<0]*sdy2,fhat2[ygrid2<0],col='green')
  lines(ygrid2[ygrid2<0]*sdy2,fhat2.u[ygrid2<0],lty=3,col='red',lwd=2)
  lines(ygrid2[ygrid2<0]*sdy2,fhat2.l[ygrid2<0],lty=3,col='red',lwd=2)
  
  
  


