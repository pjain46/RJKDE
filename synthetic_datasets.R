
# Synthetic Datasets

generate_datasets <- function(dataset_name, n=1000) {
  switch(dataset_name,
         ## 1. Basic Gaussian Mixtures
         # Two well-separated gaussians
         "two_gaussians" = {
           c(rnorm(n*0.5, mean=-2, sd=0.5), rnorm(n*0.5, mean=2, sd=0.5))
         },

         # Three overlapping gaussians
         "three_overlapping_gaussians" = {
           c(rnorm(n*0.3, mean=0, sd=1), rnorm(n*0.4, mean=3, sd=1.2), rnorm(n*0.3, mean=6, sd=0.8))
         },

         # 2. Skewed and Non-Gaussian Distributions
         # Log-normal + Gaussian mixture
         "skewed_mix" = {
           c(rlnorm(n*0.4, meanlog=0, sdlog=0.5),rnorm(n*0.6, mean=3, sd=0.7))
         },

         # Beta distribution mixture
         "beta_mix" = {
           c(rbeta(n*0.5, 2, 5)*10, rbeta(n*0.5, 5, 2)*10 - 5)
         },

         # 3. Multimodal Distributions with Varying Complexity
         # Close modes with different variances
         "vary_complexity_mix" = {
         c(rnorm(n*0.4, mean=0, sd=0.3), rnorm(n*0.3, mean=1, sd=0.1), rnorm(n*0.3, mean=1.5, sd=0.5))
         },

         # 4. Heavy-Tailed and Outlier-Prone Distributions
         # t-distribution + Gaussian mixture
         "t_normal_mix" = {
           c(rt(n*0.3, df=3)*2 + 1, rnorm(n*0.7, mean=4, sd=0.5))
         },

         # Cauchy distribution mixture
         "cauchy_mix" = {
           c(rcauchy(n*0.3, location=-2, scale=0.5), rcauchy(n*0.7, location=2, scale=1))
         },

         # 5. Uneven Cluster Sizes and Weights
         # Exponentially decreasing cluster sizes
         "uneven_mix" = {
           cluster_sizes <- round(n*exp(-(1:5)))
           unlist(mapply(rnorm, n=cluster_sizes, mean=seq(-4,4,length.out=5), sd=seq(0.2,0.6,length.out=5)))
         },

         # 6. Edge Cases and Stress Tests
         # Single component
         "single_normal" = {
           rnorm(n, mean=0, sd=1)
         },

         # Nearly identical components
         "near_identical_mix" = {
           c(rnorm(n*0.5, mean=0, sd=0.5),  rnorm(n*0.5, mean=0.1, sd=0.52))
         },

         # Extreme outliers
         "outliers" = {
           c(rnorm(n*0.9, mean=0, sd=1), rnorm(n*0.1, mean=10, sd=0.1))
         },

         stop("Unknown dataset name")
  )
}
