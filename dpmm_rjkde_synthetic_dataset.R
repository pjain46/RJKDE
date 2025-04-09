
# Time and Accuracy comparison of DPMM and RJKDE

library(RJKDE)
source("~/ASU Dropbox/Palak Jain/Palak_Hahn/Chapter2/dpmm_function.R")
source("~/ASU Dropbox/Palak Jain/Palak_Hahn/Chapter2/metrics.R")

# set.seed(123)

mc <- 20
gridsize <- 50

# Estimated Density
fhat_rjkde <- matrix(0, gridsize, mc)
fhat_dpmm <- matrix(0, gridsize, mc)

fhat_rjkde.u <- matrix(0, gridsize, mc)
fhat_dpmm.u <- matrix(0, gridsize, mc)

fhat_rjkde.l <- matrix(0, gridsize, mc)
fhat_dpmm.l <- matrix(0, gridsize, mc)

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
cc_rjkde <- rep(0, mc)
cc_dpmm <- rep(0, mc)

# Credible Interval Length
cil_rjkde <- rep(0, mc)
cil_dpmm <- rep(0, mc)




# Generate data
n <- 5000
kprime <- 4
mu <- 2 * rbeta(kprime, 0.8, 0.8) - 0.5
sig <- seq(0.125, 0.75, length.out = kprime)
ind <- sample(1:kprime, n, replace = TRUE)
x <- mu[ind] + sig[ind] * rnorm(n)
data <- sample(x, n)
ygrid <- seq(min(data) - 1, max(data) + 1, length.out = 100)



# # Generate data
# true_means <- c(-3, 0, 3)  # True cluster means
# true_sds <- c(1, 0.5, 1)    # True cluster standard deviations
# n_points <- c(50, 100, 50)  # Points per cluster
# n <- sum(n_points)
# data <- unlist(mapply(rnorm, n_points, true_means, true_sds))
# data <- data[sample(length(data))]  # Shuffling of Data points
# ygrid <- seq(min(data) - 1, max(data) + 1, length.out = 100)


# Calculate true weights (proportion in each component)
true_weights <- table(ind)/n
true_weights <- as.numeric(true_weights)  # Convert to numeric vector

# Calculate true density at each grid point
true_density <- numeric(length(ygrid))

for (i in 1:length(ygrid)) {
  # Sum weighted densities from all components
  true_density[i] <- sum(true_weights * dnorm(ygrid[i], mean = mu, sd = sig))
}

n_iter <- 100

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

# Plot results
hist(data, freq = FALSE, ylim = c(0, 1), breaks = 35, main = '', ylab = 'data')
lines(ygrid,fhat_rjkde,col='green',lwd=2)
lines(ygrid,fhat_rjkde.u,lty=3,lwd=2, col='red')
lines(ygrid,fhat_rjkde.l,lty=3,lwd=2, col='red')
lines(ygrid,fhat_dpmm,col='blue',lwd=2)
lines(ygrid,fhat_dpmm.u,lty=3,lwd=2, col='pink')
lines(ygrid,fhat_dpmm.l,lty=3,lwd=2, col='pink')


# Comparison Matrix

metrics_rjkde <- metrics(true_density, fhat_rjkde, fhat_rjkde.l, fhat_rjkde.u)
metrics_dpmm <- metrics(true_density, fhat_dpmm, fhat_dpmm.l, fhat_dpmm.u)
results_df <- data.frame(
  Mean_Time = c(elapsed_time_rjkde, elapsed_time_dpmm),
  Mean_RMSE = c(metrics_rjkde$Mean_RMSE, metrics_dpmm$Mean_RMSE),
  Mean_Coverage = c(metrics_rjkde$Mean_Coverage,metrics_dpmm$Mean_Coverage),
  Mean_Complete_Coverage = c(metrics_rjkde$Mean_Complete_Coverage,metrics_dpmm$Mean_Complete_Coverage),
  Mean_Interval_Length = c(metrics_rjkde$Mean_Interval_Length, metrics_dpmm$Mean_Interval_Length)

)
row.names(results_df)<- c("RJKDE", "DPMM")

# Print the data frame
print(round(results_df,4))
