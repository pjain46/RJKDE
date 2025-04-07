

# Metrics

metrics <- function(f_true, f_calculated, f_calculated.l, f_calculated.u, n_iter = 500) {
  # Args:
  #   f_true: True density curve on ygrid
  #   f_calculated: Estimated density on ygrid by dpmm or rjkde
  #   f_calculated.l: Lower limit of credible interval
  #   f_calculated.u: Upper limit of credible interval
  #   n_iter: Number of Gibbs sampling iterations
  #
  # Returns:
  #   A list containing:
  #     - Mean_RMSE: Root Mean Square Average on mc simulations
  #     - Mean_Coverage: How much density lie between the credible interval
  #     - Mean_Complete_Coverage: if complete density curve lies between the credible interval
  #     - Mean_Interval_Length: Credible interval width
  #

  fmax <- max(f_true)

  # Error
  err <- rep(0, n_iter)

  # Coverage
  coverage <- rep(0, n_iter)

  # Complete Coverage
  comp_coverage <- rep(0, n_iter)

  # Credible Interval Length
  credible_interval_len <- rep(0, n_iter)

  for (i in 1:n_iter){
    # Error
    err[iter] <- sqrt(mean((f_true - f_calculated)^2))/fmax

    # Coverage
    coverage[iter] <- mean(f_calculated.u > f_true & f_calculated.l <= f_true)

    # Complete Coverage
    comp_coverage[iter] <- mean(all(f_calculated.u > f_true & f_calculated.l <= f_true))

    # Credible Interval Length
    comp_coverage[iter] <- mean(f_calculated.u - f_calculated.l)

  }

  # Return results
  list(
    Mean_RMSE = mean(err),
    Mean_Coverage = mean(coverage),
    Mean_Complete_Coverage = mean(comp_coverage),
    Mean_Interval_Length = mean(comp_coverage)
  )

}
