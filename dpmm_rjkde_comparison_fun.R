
library(RJKDE)

compare_methods <- function(dataset_name, n_iter=100, mc=10, gridsize=100) {

  # Generate Data
  data <- generate_datasets(dataset_name)
  ygrid <- seq(min(data)-1, max(data)+1, length.out=gridsize)

  # True Density
  true_density <- density(data, n = length(ygrid), from = min(data) - 1, to = max(data) + 1)$y


  # Initialize storage of results
  results <- list(
    fhat_rjkde_mc = matrix(0, gridsize, mc),
    fhat_dpmm_mc = matrix(0, gridsize, mc),
    fhat_rjkde_mc.u = matrix(0, gridsize, mc),
    fhat_dpmm_mc.u = matrix(0, gridsize, mc),
    fhat_rjkde_mc.l = matrix(0, gridsize, mc),
    fhat_dpmm_mc.l = matrix(0, gridsize, mc),
    time = data.frame(RJKDE=rep(0,mc), DPMM=rep(0,mc)),
    metrics = list()
  )

  for (iter in 1:mc){

    flag <- iter%%2==0
    if(flag) print(iter)

    # Run DPMM
    # Start timer
    start_time <- Sys.time()
    dpmm <- dpmm_clustering(data, ygrid, n_iter = n_iter)
    # End timer
    end_time <- Sys.time()
    elapsed_time_dpmm <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Run RJKDE
    # Start timer
    start_time <- Sys.time()
    rjkde <- rj_mcmc_rcpp(data, ygrid, mc = n_iter)
    # End timer
    end_time <- Sys.time()
    elapsed_time_rjkde <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Store results
    results$fhat_rjkde_mc[,iter] <- rowMeans(rjkde$fsamps)
    results$fhat_dpmm_mc[,iter] <- rowMeans(dpmm$density_history)

    results$fhat_rjkde_mc.u[,iter] <- apply(rjkde$fsamps, 1, quantile, 0.95)
    results$fhat_dpmm_mc.u[,iter] <- apply(dpmm$density_history, 1, quantile, 0.95)

    results$fhat_rjkde_mc.l[,iter] <- apply(rjkde$fsamps, 1, quantile, 0.05)
    results$fhat_dpmm_mc.l[,iter] <- apply(dpmm$density_history, 1, quantile, 0.05)

    results$time$RJKDE[iter] <- elapsed_time_rjkde
    results$time$DPMM[iter] <- elapsed_time_dpmm

    # Calculate metrics
    results$metrics$rjkde[[iter]] <- metrics(true_density,
                                             results$fhat_rjkde_mc[,iter],
                                             results$fhat_rjkde_mc.l[,iter],
                                             results$fhat_rjkde_mc.u[,iter])

    results$metrics$dpmm[[iter]] <- metrics(true_density,
                                            results$fhat_dpmm_mc[,iter],
                                            results$fhat_dpmm_mc.l[,iter],
                                            results$fhat_dpmm_mc.u[,iter])

  }

  # Plot results
  hist(data, freq = FALSE, breaks = 35, main = '', ylab = 'density', xlab = 'data')
  lines(true_density,col='black',lwd=2)
  lines(ygrid,results$fhat_rjkde_mc[,iter],col='green',lwd=2)
  polygon(c(ygrid, rev(ygrid)),
          c(results$fhat_rjkde_mc.u[,iter], rev(results$fhat_rjkde_mc.l[,iter])),
          col=rgb(0,1,0,0.2), border=NA)
  # lines(ygrid,results$fhat_rjkde_mc.u[,iter],lty=3,lwd=2, col='red')
  # lines(ygrid,results$fhat_rjkde_mc.l[,iter],lty=3,lwd=2, col='red')
  lines(ygrid,results$fhat_dpmm_mc[,iter],col='blue',lwd=2)
  polygon(c(ygrid, rev(ygrid)),
          c(results$fhat_dpmm_mc.u[,iter], rev(results$fhat_dpmm_mc.l[,iter])),
          col=rgb(0,1,0,0.2), border=NA)
  # lines(ygrid,results$fhat_dpmm_mc.u[,iter],lty=3,lwd=2, col='orange')
  # lines(ygrid,results$fhat_dpmm_mc.l[,iter],lty=3,lwd=2, col='orange')
  legend("topright", legend=c("Truth", "RJKDE", "DPMM"), col=c("black", "green", "blue"), lwd=2, , cex = 0.6)


  return(results)
}


generate_results_table <- function(results) {
  # Calculate mean metrics across all iterations
  mean_metrics <- data.frame(
    Method = c("RJKDE", "DPMM"),
    Time = c(mean(results$time$RJKDE), mean(results$time$DPMM)),
    RMSE = c(mean(sapply(results$metrics$rjkde, function(x) x$Mean_RMSE)),
             mean(sapply(results$metrics$dpmm, function(x) x$Mean_RMSE))),
    Coverage = c(mean(sapply(results$metrics$rjkde, function(x) x$Mean_Coverage)),
                 mean(sapply(results$metrics$dpmm, function(x) x$Mean_Coverage))),
    Complete_Coverage = c(mean(sapply(results$metrics$rjkde, function(x) x$Mean_Complete_Coverage)),
                          mean(sapply(results$metrics$dpmm, function(x) x$Mean_Complete_Coverage))),
    Interval_Width = c(mean(sapply(results$metrics$rjkde, function(x) x$Mean_Interval_Width)),
                       mean(sapply(results$metrics$dpmm, function(x) x$Mean_Interval_Width)))
  )

  # Round only numeric columns
  numeric_cols <- sapply(mean_metrics, is.numeric)
  mean_metrics[numeric_cols] <- round(mean_metrics[numeric_cols], 4)

  return(mean_metrics)
}


