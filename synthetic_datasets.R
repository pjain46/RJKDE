
# Synthetic Datasets

generate_datasets <- function(dataset_name, n=1000) {
  switch(dataset_name,
         ## 1. Basic Gaussian Mixtures
         # Two well-separated gaussians
         "two_gaussians" = {
           data <- c(rnorm(n * 0.5, mean = -2, sd = 0.5),
                     rnorm(n * 0.5, mean = 2, sd = 0.5))
           density <- function(x) {
             0.5 * dnorm(x, mean = -2, sd = 0.5) +
               0.5 * dnorm(x, mean = 2, sd = 0.5)
           }
           list(data = data, true_density = density)
         },

         # Three overlapping gaussians
         "three_overlapping_gaussians" = {
           data <- c(rnorm(n * 0.3, mean = 0, sd = 1),
                     rnorm(n * 0.4, mean = 3, sd = 1.2),
                     rnorm(n * 0.3, mean = 6, sd = 0.8))
           density <- function(x) {
             0.3 * dnorm(x, 0, 1) +
               0.4 * dnorm(x, 3, 1.2) +
               0.3 * dnorm(x, 6, 0.8)
           }
           list(data = data, true_density = density)
         },

         # 2. Skewed and Non-Gaussian Distributions
         # Log-normal + Gaussian mixture
         "skewed_mix" = {
           data <- c(rlnorm(n * 0.4, meanlog = 0, sdlog = 0.5),
                     rnorm(n * 0.6, mean = 3, sd = 0.7))
           density <- function(x) {
             0.4 * dlnorm(x, meanlog = 0, sdlog = 0.5) +
               0.6 * dnorm(x, mean = 3, sd = 0.7)
           }
           list(data = data, true_density = density)
         },

         # Beta distribution mixture
         "beta_mix" = {
           data <- c(rbeta(n * 0.5, 2, 5) * 10,
                     rbeta(n * 0.5, 5, 2) * 10 - 5)
           density <- function(x) {
             d1 <- dbeta(x / 10, 2, 5) / 10
             d2 <- dbeta((x + 5) / 10, 5, 2) / 10
             0.5 * d1 + 0.5 * d2
           }
           list(data = data, true_density = density)
         },

         # 3. Multimodal Distributions with Varying Complexity
         # Close modes with different variances
         "vary_complexity_mix" = {
           data <- c(rnorm(n * 0.4, 0, 0.3),
                     rnorm(n * 0.3, 1, 0.1),
                     rnorm(n * 0.3, 1.5, 0.5))
           density <- function(x) {
             0.4 * dnorm(x, 0, 0.3) +
               0.3 * dnorm(x, 1, 0.1) +
               0.3 * dnorm(x, 1.5, 0.5)
           }
           list(data = data, true_density = density)
         },

         # 4. Heavy-Tailed and Outlier-Prone Distributions
         # t-distribution + Gaussian mixture
         "t_normal_mix" = {
           data <- c(rt(n * 0.3, df = 3) * 2 + 1,
                     rnorm(n * 0.7, 4, 0.5))
           density <- function(x) {
             0.3 * dt((x - 1) / 2, df = 3) / 2 +
               0.7 * dnorm(x, 4, 0.5)
           }
           list(data = data, true_density = density)
         },

         # Cauchy distribution mixture
         "cauchy_mix" = {
           data <- c(rcauchy(n * 0.3, location = -2, scale = 0.5),
                     rcauchy(n * 0.7, location = 2, scale = 1))
           density <- function(x) {
             0.3 * dcauchy(x, location = -2, scale = 0.5) +
               0.7 * dcauchy(x, location = 2, scale = 1)
           }
           list(data = data, true_density = density)
         },

         # 5. Uneven Cluster Sizes and Weights
         # Exponentially decreasing cluster sizes
         "uneven_mix" = {
           cluster_sizes <- round(n * exp(-(1:5)))
           mus <- seq(-4, 4, length.out = 5)
           sds <- seq(0.2, 0.6, length.out = 5)
           total <- sum(cluster_sizes)
           data <- unlist(mapply(rnorm, n = cluster_sizes, mean = mus, sd = sds))
           weights <- cluster_sizes / total
           density <- function(x) {
             rowSums(sapply(1:5, function(i) weights[i] * dnorm(x, mean = mus[i], sd = sds[i])))
           }
           list(data = data, true_density = density)
         },

         # 6. Edge Cases and Stress Tests
         # Single component
         "single_normal" = {
           data <- rnorm(n, 0, 1)
           density <- function(x) dnorm(x, 0, 1)
           list(data = data, true_density = density)
         },

         # Nearly identical components
         "near_identical_mix" = {
           data <- c(rnorm(n * 0.5, 0, 0.5),
                     rnorm(n * 0.5, 0.1, 0.52))
           density <- function(x) {
             0.5 * dnorm(x, 0, 0.5) +
               0.5 * dnorm(x, 0.1, 0.52)
           }
           list(data = data, true_density = density)
         },

         # Extreme outliers
         "outliers" = {
           data <- c(rnorm(n * 0.9, 0, 1),
                     rnorm(n * 0.1, 10, 0.1))
           density <- function(x) {
             0.9 * dnorm(x, 0, 1) +
               0.1 * dnorm(x, 10, 0.1)
           }
           list(data = data, true_density = density)
         },

         stop("Unknown dataset name")
  )
}
