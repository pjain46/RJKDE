
# Gibbs Sampler for Dirichlet Process Mixture Model (DPMM)
# Assumes Gaussian likelihood with unknown mean/variance

library(ggplot2)

set.seed(123)

# Generate Synthetic Data
true_means <- c(-3, 0, 3)  # True cluster means
true_sds <- c(1, 0.5, 1)    # True cluster standard deviations
n_points <- c(50, 100, 50)  # Points per cluster

data <- unlist(mapply(rnorm, n_points, true_means, true_sds))
data <- data[sample(length(data))]  # Shuffling of Data points


### 2. Initialize Gibbs Sampler ###
n_iter <- 500
n <- length(data)

# Drichlet Process Concentration Parameter
alpha <- 1

# Hyper parameters
tau <- 1 # Known Precision for Normal - Normal conjugacy model

mu0 <- 0 # Prior mean for cluster means
tau0 <- 0.1 # Precision of prior for means
a0 <- 2 # Shape for precision prior
b0 <- 1 # Scale for precision prior


# All data points in one cluster -- Initialize cluster assignment
z = rep(1, n)
k = 1 # Number of clusters

# Store results
cluster_history <- matrix(0, n_iter, n)
n_clusters <- numeric(n_iter)

for (iter in 1:n_iter){

  # Step 1: Update cluster assignments -- Chinese restaurant process

  for (i in 1:k){

    # If more than one point is in the cluster
    # Temporarily unassign point i but keep the cluster alive

    if(sum(z == z[i]) > 1){
      z[i] <- NA
    }else{

      # If point i is the ONLY member of its cluster
      # Decrease total cluster count since we're deleting this cluster
      # Renumber higher clusters to fill the gap
      # Unassign point i
      k <- k-1
      z[z > z[i]] <- z[z > z[i]] - 1
      z[i] <- NA
    }

    # Compute probabilities for existing clusters -- using Normal - Normal conjugacy model
    cluster_probs <- sapply(1:k, function(j){
      nj = sum(z == j, na.rm = TRUE)
      x_bar <- mean(data[z ==j], na.rm = TRUE)
      prec_j <- tau0 + nj*tau
      mu_j <- (tau0*mu0 + tau*nj*x_bar)/prec_j

      # Compute probability for cluster j: Likelihood × Cluster size
      # Posterior Predictive Distribution from conjugate model
      dnorm(data[i], mu_j, sqrt(1/prec_j + 1/tau)) * nj
    })

    new_cluster_prob <- alpha * dnorm(data[i], mu0, sqrt(1/tau0 + 1/tau))


    # Sample new assignment for the data point
    probs <- c(cluster_probs, new_cluster_prob)
    z[i] <- sample(1:(k+1), 1, prob = probs)

    # If we sampled the new cluster
    if (z[i] == k+1){
      k  <- k+1
    }
  }

  # Step 2: Update cluster parameters
  cluster_means <- rep(0, k)
  cluster_precisions <- rep(0,k)

  for (j in 1:k){
    nj <- sum(z == j)
    x_bar <- mean(data[z == j])

    # Posterior parameters

    # Posterior precision = Prior precision + Data precision
    prec_j <- tau0 + nj*tau
    mu_j <- (tau0*mu0 + tau*nj*x_bar)/prec_j

    a_j <- a0 + nj/2
    # b_new = b_0 + 0.5 * SSE within cluster + 0.5 * Penalty from deviating from mu_0
    b_j <- b0 + 0.5 * sum((data[z == j] - x_bar)^2) +
      (tau0 * nj * (x_bar - mu0)^2) / (2 * (tau0 + nj))

    # Sample new parameters
    cluster_precisions[j] <- rgamma(1, a_j, b_j)
    cluster_means[j] <- rnorm(1, mu_j, sqrt(1/(prec_j * cluster_precisions[j])))
  }

  # Store results
  cluster_history[iter, ] <- z
  n_clusters[iter] <- k

}




