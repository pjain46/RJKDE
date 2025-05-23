

# Metrics

metrics <- function(f_true, f_est, f_lower, f_upper) {
  # Args:
  #   f_true: True density curve on ygrid
  #   f_est: Estimated density on ygrid by dpmm or rjkde
  #   f_lower: Lower limit of credible interval
  #   f_upper: Upper limit of credible interval
  #
  # Returns:
  #   A list containing:
  #     - Mean_RMSE: Root Mean Square Average on mc simulations
  #     - Mean_Coverage: How much density lie between the credible interval
  #     - Mean_Complete_Coverage: if complete density curve lies between the credible interval
  #     - Mean_Interval_Length: Credible interval width
  #

  fmax <- max(f_true)
  if (fmax == 0) fmax <- 1e-8

  # Normalized RMSE
  rmse <- sqrt(mean((f_true - f_est)^2)) / fmax

  # Pointwise coverage
  coverage <- mean(f_upper >= f_true & f_lower <= f_true)

  # Complete coverage (entire curve within CI)
  complete_coverage <- as.numeric(all(f_upper >= f_true & f_lower <= f_true))

  # Average interval width
  interval_width <- mean(f_upper - f_lower)

  # Return results
  list(
    Mean_RMSE = rmse,
    Mean_Coverage = coverage,
    Mean_Complete_Coverage = complete_coverage,
    Mean_Interval_Width = interval_width
  )

}
