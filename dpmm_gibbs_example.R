
# Gibbs Sampler for Dirichlet Process Mixture Model (DPMM)
# Assumes Gaussian likelihood with unknown mean/variance

# For Dirichlet and Inverse-Gamma distributions: MCMCpack library
library(MCMCpack)
library(ggplot2)

set.seed(123)

# Generate Synthetic Data
true_means <- c(-3, 0, 3)  # True cluster means
true_sds <- c(1, 0.5, 1)    # True cluster standard deviations
n_points <- c(50, 100, 50)  # Points per cluster

data <- unlist(mapply(rnorm, n_points, true_means, true_sds))
data <- data[sample(length(data))]  # Shuffling of Data points


