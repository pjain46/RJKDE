library(Rcpp)


# generate data

n <- 5000

kprime = 4
mu = 2*rbeta(kprime, 0.8,0.8)-0.5
sig = seq(0.125,0.75, length.out = kprime)

ind <- sample(1:kprime,n,replace=TRUE)
x = mu[ind] + sig[ind]*rnorm(n)
data <- sample(x,n,replace=TRUE)

hist(data,freq=FALSE,ylim=c(0,1), 35,main='',ylab='data')

### 2. Initialize Gibbs Sampler ###
n_iter <- 500
n <- length(data)

# Drichlet Process Concentration Parameter
alpha <- 1

# Hyper parameters
lambda <- 1 # Known Precision for Normal - Normal conjugacy model

mu0 <- 0 # Prior mean for cluster means
lambda0 <- 0.1 # Precision of prior for means
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

  for (i in 1:n){

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
      prec_j <- lambda0 + nj*lambda
      mu_j <- (lambda0*mu0 + lambda*nj*x_bar)/prec_j

      # Compute probability for cluster j: Likelihood × Cluster size
      # Posterior Predictive Distribution from conjugate model
      dnorm(data[i], mu_j, sqrt(1/prec_j + 1/lambda)) * nj
    })

    new_cluster_prob <- alpha * dnorm(data[i], mu0, sqrt(1/lambda0 + 1/lambda))


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
    prec_j <- lambda0 + nj*lambda
    mu_j <- (lambda0*mu0 + lambda*nj*x_bar)/prec_j

    a_j <- a0 + nj/2
    # b_new = b_0 + 0.5 * SSE within cluster + 0.5 * Penalty from deviating from mu_0
    b_j <- b0 + 0.5 * sum((data[z == j] - x_bar)^2) +
      (lambda0 * nj * (x_bar - mu0)^2) / (2 * (lambda0 + nj))

    # Sample new parameters
    cluster_precisions[j] <- rgamma(1, a_j, b_j)
    cluster_means[j] <- rnorm(1, mu_j, sqrt(1/(prec_j * cluster_precisions[j])))
  }

  # Store results
  cluster_history[iter, ] <- z
  n_clusters[iter] <- k

}


# Plot number of clusters over iterations
plot(1:n_iter, n_clusters, type = "l",
     xlab = "Iteration", ylab = "Number of clusters",
     main = "DPMM Gibbs Sampler: Cluster Count")

# Posterior cluster assignments (last iteration)
df <- data.frame(x = data, cluster = factor(z))
ggplot(df, aes(x, fill = cluster)) +
  geom_histogram(binwidth = 0.5, alpha = 0.7) +
  ggtitle("DPMM Clustering Results") +
  theme_minimal()
