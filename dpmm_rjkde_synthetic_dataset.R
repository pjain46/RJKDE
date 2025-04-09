
# Time and Accuracy comparison of DPMM and RJKDE

library(RJKDE)
source("~/ASU Dropbox/Palak Jain/Palak_Hahn/Chapter2/dpmm_function.R")
source("~/ASU Dropbox/Palak Jain/Palak_Hahn/Chapter2/metrics.R")

# set.seed(123)

# Number of iterations for DPMM and RJKDE
n_iter <- 100

#Monte Carlo Setup:
mc <- 10
gridsize <- 100

# Estimated Density
fhat_rjkde_mc <- matrix(0, gridsize, mc)
fhat_dpmm_mc <- matrix(0, gridsize, mc)

fhat_rjkde_mc.u <- matrix(0, gridsize, mc)
fhat_dpmm_mc.u <- matrix(0, gridsize, mc)

fhat_rjkde_mc.l <- matrix(0, gridsize, mc)
fhat_dpmm_mc.l <- matrix(0, gridsize, mc)

# Time
t_rjkde <- rep(0, mc)
t_dpmm <- rep(0, mc)

# Error
err_rjkde <- rep(0, mc)
err_dpmm <- rep(0, mc)

# Coverage
cover_rjkde<- rep(0, mc)
cover_dpmm <- rep(0, mc)

# Complete Coverage
complete_coverage_rjkde <- rep(0, mc)
complete_coverage_dpmm <- rep(0, mc)

# Credible Interval Length
interval_width_rjkde <- rep(0, mc)
interval_width_dpmm <- rep(0, mc)



for (iter in 1:mc){

  flag <- iter%%2==0
  if(flag) print(iter)

  # Generate data
  # n <- 5000
  # kprime <- 4
  # mu <- 2 * rbeta(kprime, 0.8, 0.8) - 0.5
  # sig <- seq(0.125, 0.75, length.out = kprime)
  # ind <- sample(1:kprime, n, replace = TRUE)
  # x <- mu[ind] + sig[ind] * rnorm(n)
  # data <- sample(x, n)
  # ygrid <- seq(min(data) - 1, max(data) + 1, length.out = gridsize)



  # # Generate data
  # true_means <- c(-3, 0, 3)  # True cluster means
  # true_sds <- c(1, 0.5, 1)    # True cluster standard deviations
  # n_points <- c(50, 100, 50)  # Points per cluster
  # n <- sum(n_points)
  # data <- unlist(mapply(rnorm, n_points, true_means, true_sds))
  # data <- data[sample(length(data))]  # Shuffling of Data points
  # ygrid <- seq(min(data) - 1, max(data) + 1, length.out = gridsize)

  # Two well-separated gaussians
  data <- c(rnorm(500, mean=-2, sd=0.5), rnorm(500, mean=2, sd=0.5))
  ygrid <- seq(min(data) - 1, max(data) + 1, length.out = gridsize)

  #  t-distribution + Gaussian mixture
  data <- c(rt(300, df=3)*2 + 1, rnorm(700, mean=4, sd=0.5))
  ygrid <- seq(min(data) - 1, max(data) + 1, length.out = gridsize)
  n <- length(data)


  # # Calculate true weights (proportion in each component)
  # true_weights <- table(ind)/n
  # true_weights <- as.numeric(true_weights)  # Convert to numeric vector
  #
  # # Calculate true density at each grid point
  # true_density <- numeric(length(ygrid))
  #
  # for (i in 1:length(ygrid)) {
  #   # Sum weighted densities from all components
  #   true_density[i] <- sum(true_weights * dnorm(ygrid[i], mean = mu, sd = sig))
  # }



  # Run DPMM clustering
  # Start timer
  start_time <- Sys.time()
  results <- dpmm_clustering(data, ygrid, n_iter = n_iter)
  # End timer
  end_time <- Sys.time()
  elapsed_time_dpmm <- as.numeric(difftime(end_time, start_time, units = "secs"))

  fhat_dpmm <- rowMeans(results$density_history)
  fhat_dpmm.u <- apply(results$density_history,1,quantile,0.95)
  fhat_dpmm.l <- apply(results$density_history,1,quantile,0.05)

  # Hyperparameters for RJKDE
  # kinit <- 100
  # burn <- 1
  # kde <- density(sample(data,kinit))
  # bw <- kde$bw


  # Run RJKDE
  # Start timer
  start_time <- Sys.time()
  KDEpost <- rj_mcmc_rcpp(data, ygrid, mc = n_iter)
  # End timer
  end_time <- Sys.time()
  elapsed_time_rjkde <- as.numeric(difftime(end_time, start_time, units = "secs"))

  fsamps <- KDEpost$fsamps
  fhat_rjkde <- rowMeans(fsamps)
  fhat_rjkde.u <- apply(fsamps,1,quantile,0.95)
  fhat_rjkde.l <- apply(fsamps,1,quantile,0.05)


  # Estimated Density
  fhat_rjkde_mc[,iter] <- fhat_rjkde
  fhat_dpmm_mc[,iter] <- fhat_dpmm

  fhat_rjkde_mc.u[,iter] <- fhat_rjkde.u
  fhat_dpmm_mc.u[,iter] <- fhat_dpmm.u

  fhat_rjkde_mc.l[,iter] <- fhat_rjkde.l
  fhat_dpmm_mc.l[,iter] <- fhat_dpmm.l

  # Comparison Metrics
  metrics_rjkde <- metrics(true_density, fhat_rjkde, fhat_rjkde.l, fhat_rjkde.u)
  metrics_dpmm <- metrics(true_density, fhat_dpmm, fhat_dpmm.l, fhat_dpmm.u)


  # Time
  t_rjkde[iter] <- elapsed_time_rjkde
  t_dpmm[iter] <- elapsed_time_dpmm

  # Error
  err_rjkde[iter] <- metrics_rjkde$Mean_RMSE
  err_dpmm[iter] <- metrics_dpmm$Mean_RMSE

  # Coverage
  cover_rjkde[iter] <- metrics_rjkde$Mean_Coverage
  cover_dpmm[iter] <- metrics_dpmm$Mean_Coverage

  # Complete Coverage
  complete_coverage_rjkde[iter] <- metrics_rjkde$Mean_Complete_Coverage
  complete_coverage_dpmm[iter] <- metrics_dpmm$Mean_Complete_Coverage

  # Credible Interval Length
  interval_width_rjkde[iter] <- metrics_rjkde$Mean_Interval_Width
  interval_width_dpmm[iter] <- metrics_dpmm$Mean_Interval_Width

}

# Plot results
hist(data, freq = FALSE, ylim = c(0, 1), breaks = 35, main = '', ylab = 'density', xlab = 'data')
lines(density(data),col='black',lwd=2)
lines(ygrid,fhat_rjkde,col='green',lwd=2)
lines(ygrid,fhat_rjkde.u,lty=3,lwd=2, col='red')
lines(ygrid,fhat_rjkde.l,lty=3,lwd=2, col='red')
lines(ygrid,fhat_dpmm,col='blue',lwd=2)
lines(ygrid,fhat_dpmm.u,lty=3,lwd=2, col='orange')
lines(ygrid,fhat_dpmm.l,lty=3,lwd=2, col='orange')
legend("topright", legend=c("Truth", "RJKDE", "DPMM"), col=c("black", "green", "blue"), lwd=2, , cex = 0.6)


# Comparison Metrices
results_df <- data.frame(
  Mean_Time = c(mean(t_rjkde), mean(t_dpmm)),
  Mean_RMSE = c(mean(err_rjkde), mean(err_dpmm)),
  Mean_Coverage = c(mean(cover_rjkde),mean(cover_dpmm)),
  Mean_Complete_Coverage = c(mean(complete_coverage_rjkde),mean(complete_coverage_dpmm)),
  Mean_Interval_Width = c(mean(interval_width_rjkde), mean(interval_width_dpmm))

)
row.names(results_df)<- c("RJKDE", "DPMM")

# Print the data frame
print(round(results_df,4))
