# run_analysis.R
source("~/ASU Dropbox/Palak Jain/Palak_Hahn/Chapter2/synthetic_datasets.R")
source("~/ASU Dropbox/Palak Jain/Palak_Hahn/Chapter2/dpmm_rjkde_comparison_fun.R")
source("~/ASU Dropbox/Palak Jain/Palak_Hahn/Chapter2/dpmm_function.R")
source("~/ASU Dropbox/Palak Jain/Palak_Hahn/Chapter2/metrics.R")

# List of datasets to analyze
datasets <- c("two_gaussians", "three_overlapping_gaussians", "skewed_mix",
              "beta_mix", "vary_complexity_mix", "t_normal_mix", "cauchy_mix", "uneven_mix",
              "single_normal", "near_identical_mix", "outliers")

# Run analysis
dataset <- datasets[1]
message("\nAnalyzing dataset: ", dataset)
results <- compare_methods(dataset, n_iter=100, mc=5)

# Generate and save results table
results_table <- generate_results_table(results)
print(results_table)

