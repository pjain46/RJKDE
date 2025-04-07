
# Dirichlet Process Mixture Model (DPMM) clustering using Gibbs sampling


dpmm_clustering <- function(data, ygrid = seq(-3, 3, length.out = 100),
                            n_iter = 500,
                            alpha = 1, lambda = 1,
                            mu0 = 0, lambda0 = 0.1,
                            a0 = 2, b0 = 1) {

  # Args:
  #   data: Numeric vector of data to cluster
  #   ygrid: Grid points where density should be evaluated
  #   n_iter: Number of Gibbs sampling iterations
  #   alpha: DP concentration parameter
  #   lambda: Known precision for Normal-Normal conjugacy model
  #   mu0: Prior mean for cluster means for Normal-Normal conjugacy model
  #   lambda0: Precision of prior for means
  #   a0: Shape for precision prior for Normal - (Normal - Gamma) conjugacy model for mean and precision to be assumed exchangeable
  #   b0: Scale for precision prior
  #
  # Returns:
  #   A list containing:
  #     - cluster_history: Matrix of cluster assignments for each iteration
  #     - n_clusters: Vector of number of clusters at each iteration
  #     - parameters: List of final cluster means and precisions


  n <- length(data)
  ngrid <- length(ygrid)

  # Initialize cluster assignments - all data points in one cluster
  z <- rep(1, n)
  k <- 1  # Number of clusters

  # Store results
  cluster_history <- matrix(0, n_iter, n)
  n_clusters <- numeric(n_iter)
  density_history <- matrix(0, ngrid, n_iter)

  for (iter in 1:n_iter) {
    # Step 1: Update cluster assignments -- Chinese restaurant process
    for (i in 1:n) {

      # If more than one data point is in the cluster
      # Temporarily unassign point i
      if (sum(z == z[i]) > 1) {
        z[i] <- NA
      } else {
        # If point i is the ONLY member of its cluster
        # Decrease total cluster count since we're deleting this cluster
        # Renumber higher clusters to fill the gap
        # Unassign point i
        k <- k - 1
        z[z > z[i]] <- z[z > z[i]] - 1
        z[i] <- NA
      }

      # Compute probabilities for existing clusters -- using Normal - Normal conjugacy model
      cluster_probs <- sapply(1:k, function(j) {
        nj <- sum(z == j, na.rm = TRUE)
        x_bar <- mean(data[z == j], na.rm = TRUE)
        prec_j <- lambda0 + nj * lambda
        mu_j <- (lambda0 * mu0 + lambda * nj * x_bar) / prec_j

        # Compute probability for cluster j: Likelihood × Cluster size
        # Posterior Predictive Distribution from conjugate model
        dnorm(data[i], mu_j, sqrt(1 / prec_j + 1 / lambda)) * nj
      })

      # Probability for new cluster
      new_cluster_prob <- alpha * dnorm(data[i], mu0, sqrt(1 / lambda0 + 1 / lambda))

      # Sample new assignment for the data point
      probs <- c(cluster_probs, new_cluster_prob)
      z[i] <- sample(1:(k + 1), 1, prob = probs)

      # If we sampled the new cluster
      if (z[i] == k + 1) {
        k <- k + 1
      }
    }

    # Step 2: Update cluster parameters
    cluster_means <- numeric(k)
    cluster_precisions <- numeric(k)

    for (j in 1:k) {
      nj <- sum(z == j)
      x_bar <- mean(data[z == j])

      # Posterior parameters
      # Posterior precision = Prior precision + Data precision
      prec_j <- lambda0 + nj * lambda
      mu_j <- (lambda0 * mu0 + lambda * nj * x_bar) / prec_j

      a_j <- a0 + nj / 2
      # b_new = b_0 + 0.5 * SSE within cluster + 0.5 * Penalty from deviating from mu_0
      b_j <- b0 + 0.5 * sum((data[z == j] - x_bar)^2) +
        (lambda0 * nj * (x_bar - mu0)^2) / (2 * (lambda0 + nj))

      # Sample new parameters
      cluster_precisions[j] <- rgamma(1, a_j, b_j)
      cluster_means[j] <- rnorm(1, mu_j, sqrt(1 / (prec_j * cluster_precisions[j])))
    }

    # Store results
    cluster_history[iter, ] <- z
    n_clusters[iter] <- k

    # Calculate density on ygrid for this iteration
    density_iter <- rep(0, ngrid)

    if (k > 0) {
      # For each component in the mixture
      for (j in 1:k) {
        # Get component weight (proportion of points in this cluster)
        weight <- sum(z == j) / n

        # Calculate component density on grid
        component_density <- weight * dnorm(ygrid,
                                            mean = cluster_means[j],
                                            sd = 1/sqrt(cluster_precisions[j]))

        # Add to total density
        density_iter <- density_iter + component_density
      }
    }

    # Store density for this iteration
    density_history[, iter] <- density_iter

    }

  # Return results
  list(
    cluster_history = cluster_history,
    n_clusters = n_clusters,
    parameters = list(
      means = cluster_means,
      precisions = cluster_precisions
    ),
    density_history = density_history,
    ygrid = ygrid
  )
}

# Example usage:
# Generate data
set.seed(123)
n <- 5000
kprime <- 4
mu <- 2 * rbeta(kprime, 0.8, 0.8) - 0.5
sig <- seq(0.125, 0.75, length.out = kprime)
ind <- sample(1:kprime, n, replace = TRUE)
x <- mu[ind] + sig[ind] * rnorm(n)
data <- sample(x, n, replace = TRUE)
ygrid <- seq(min(data) - 1, max(data) + 1, length.out = 100)

# Calculate true weights (proportion in each component)
true_weights <- table(ind)/n
true_weights <- as.numeric(true_weights)  # Convert to numeric vector


# Calculate true density at each grid point
true_density <- numeric(length(ygrid))

for (i in 1:length(ygrid)) {
  # Sum weighted densities from all components
  true_density[i] <- sum(true_weights * dnorm(ygrid[i], mean = mu, sd = sig))
}


# Start timer
start_time <- Sys.time()

# Run DPMM clustering
n_iter = 500
results <- dpmm_clustering(data, ygrid, n_iter = n_iter)


fhat <- rowMeans(results$density_history)
fhat.u <- apply(results$density_history,1,quantile,0.95)
#fhat.l <- rep(0,length(fhat.u))
fhat.l <- apply(results$density_history,1,quantile,0.05)

# End timer
end_time <- Sys.time()
elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Plot results
hist(data, freq = FALSE, ylim = c(0, 1), breaks = 35, main = '', ylab = 'data')
lines(ygrid,fhat,col='green',lwd=2,lty=2)
lines(ygrid,fhat.u,lty=3,col='red')
lines(ygrid,fhat.l,lty=3,col='red')


metrics <- metrics(true_density, fhat, fhat.l, fhat.u, n_iter = n_iter)
results_df <- data.frame(
  Mean_Time = 0,
  Mean_RMSE = metrics$Mean_RMSE,
  Mean_Coverage = metrics$Mean_Coverage,
  Mean_Complete_Coverage = metrics$Mean_Complete_Coverage,
  Mean_Interval_Length = metrics$Mean_Interval_Length

)
#row.names(results_df)<- c("Targeted Sub-sampling", "RJKDE Uniform Sub-sampling", "Full Data", "Bayesm with Uniform Sub-sampling")

# Print the data frame
print(round(results_df,2))
