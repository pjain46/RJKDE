#include <Rcpp.h>
using namespace Rcpp;


//' Sample from a Gaussian Mixture Model
//'
//' Generates random samples from a mixture of normal distributions with specified means
//' and a common standard deviation.
//'
//' @param mu: A numeric vector of component means (e.g., c(-1, 0, 1) for a 3-component mixture).
//' @param sigma: The common standard deviation for all mixture components (must be > 0).
//' @return A numeric vector of length 1 containing a single random sample from the mixture distribution.
//'
//' @details
//' This function implements the following steps:
//' \enumerate{
//'   \item Randomly selects one component from the mixture with equal probability
//'   \item Generates a sample from N($\mu_i$, $\sigma^2$) where $\mu_i$ is the selected component's mean
//' }
//'
//' @examples
//' # Sample from a 2-component mixture
//' means <- c(-1, 1)  # Means at -1 and 1
//' rnormmix_rcpp(means, 0.5)  # Common SD = 0.5
//'
//' # Sample from a 3-component mixture
//' rnormmix_rcpp(c(-2, 0, 2), 1.0)
//'
//' @seealso \code{\link{rnorm}} for sampling from single normal distribution
//' @export
//'
// [[Rcpp::export]]
NumericVector rnormmix_rcpp(NumericVector mu, double sigma) {
  int k = mu.length();
  IntegerVector inds = sample(k, 1, true);

  NumericVector xsamp(1);
    double x = R::rnorm(mu[inds[0]-1], sigma);
    xsamp[0] = x;

  return(xsamp);
}

// NumericVector Cquantile(NumericVector x, NumericVector q) {
//   NumericVector y = clone(x);
//   std::sort(y.begin(), y.end());
//   return y[x.size()*(q - 0.000000001)];
// }


//' Compute Log-Density of a Gaussian Mixture Model
//'
//' Calculates the log-probability density at data point `x` for a mixture of normal distributions
//' with given means, standard deviations, and uniform weights.
//'
//' @param mu: A numeric vector of component means of length k.
//' @param sigma: A numeric vector of standard deviations for each component of length k.
//' @param w: A numeric vector of component weights (in this case it is a vector of length k of ones).
//' @param x: Scalar value where density is evaluated.
//'
//' @return The log-density value at `x`.
//'
//' @details
//' The density is computed as:
//' \deqn{\log\left(\sum_{i=1}^k \frac{w_i}{k} \cdot \phi(x|\mu_i,\sigma_i)\right)}
//' where \eqn{\phi(x|\mu_i,\sigma_i)} is the normal PDF with mean \eqn{\mu_i} and standard deviation \eqn{\sigma_i}.
//'
//' @note
//' Weights are automatically normalized by the number of components (k), so the input weights
//' need not sum to 1. For standard weighted mixtures, pass pre-normalized weights.
//'
//' @examples
//' # Evaluate at x=0 for 2-component mixture
//' KDE_density_rcpp(mu = c(-1, 1),
//'                  sigma = c(0.5, 0.5),
//'                  w = c(1, 1),  # Weights automatically normalized
//'                  x = 0)
//'
//' # Compare with single normal distribution
//' KDE_density_rcpp(mu = 0, sigma = 1, w = 1, x = 0)  # log(dnorm(0,0,1))
//'
//' @seealso \code{\link{dnorm}} for single-component density
//' @export
// [[Rcpp::export]]
double KDE_density_rcpp(NumericVector mu, NumericVector sigma, NumericVector w, double x){
  int k = mu.length();
  NumericVector A(k);
  for(int i = 0; i < k; i++){
    A[i] = ((R::dnorm(x, mu[i], sigma[i], false))*w[i])/k;
  }
  double a = log(sum(A));
  return(a);
}

//' Run Reversible Jump Markov Chain Monte Carlo for Gaussian Mixture Models
//'
//' Performs Bayesian estimation of Gaussian mixture models with variable number of components
//' using Reversible Jump MCMC (RJMCMC).
//'
//' @param y: A numeric vector of observed data points.
//' @param ygrid: A numeric vector of grid points for density estimation (eg., ygrid <- seq(min(y),max(y),length.out = 500)).
//' @param drop_add_prob Numeric vector of probabilities for (drop, keep, add) moves
//'        (default: c(0.45, 0.1, 0.45)).
//' @param sig_a: Shape parameter for Beta prior on bandwidths (default: 5).
//' @param sig_b: Scale parameter for Beta prior on bandwidths (default: 5).
//' @param mu_prior_b: Prior parameter for means (default: 3).
//' @param mu_step_size: Step size for mean proposals (default: 0.1).
//' @param k: Initial number of mixture components (default: 20).
//' @param bw: Bandwidth parameter (default: 0.2).
//' @param mc: Number of MCMC iterations (default: 5000).
//' @param mu_h: Hyperparameter for mean adjustment (default: 0.0).
//' @param sig_h: Hyperparameter for bandwidth adjustment (default: 1.0).
//'
//' @return A List object containing:
//' \describe{
//'   \item{fsamps}{Matrix of density samples (columns) evaluated at ygrid points (rows)}
//'   \item{fsamps_adj}{Adjusted density samples}
//'   \item{ksamps}{Vector of component counts across iterations}
//'   \item{ygrid}{Input grid points (for reference)}
//' }
//'
//' @details
//' The algorithm implements:
//' \enumerate{
//'   \item Component add/drop component moves via RJMCMC
//'   \item Adaptive proposals for means and bandwidths
//'   \item Posterior density estimation on provided grid
//' }
//'
//' @examples
//' \dontrun{
//' # Generate sample data
//' set.seed(123)
//' y <- c(rnorm(200, -2, 1), rnorm(300, 2, 1))
//'
//' # Create evaluation grid
//' ygrid <- seq(-5, 5, length.out = 100)
//'
//' # Run MCMC
//' results <- rj_mcmc_rcpp(y, ygrid, mc = 1000)
//'
//' # Plot posterior mean density
//' plot(ygrid, rowMeans(results$fsamps), type = "l",
//'      main = "Posterior Density Estimate")
//' }
//'
//' @seealso \code{\link{KDE_density_rcpp}} for the density calculation function
//' @export
// [[Rcpp::export]]
List rj_mcmc_rcpp(NumericVector y, NumericVector ygrid, NumericVector drop_add_prob = NumericVector::create(0.45,0.1,0.45), double sig_a = 5, double sig_b = 5, double mu_prior_b = 3,double mu_step_size = 0.1, int k = 20, double bw = 0.2, int mc = 5000, double mu_h = 0.0, double sig_h = 1.0){

  double nsamp = y.length();
  double gridsize = ygrid.length();

  NumericMatrix fsamps(gridsize, mc);


  NumericMatrix fsamps_adj(gridsize, mc);


  NumericVector ksamps(mc);

  int mu_length = 2*nsamp + 500;
  NumericVector mu(mu_length);


   y.sort();

   NumericVector w(mu_length,1.0);

  int step = floor(nsamp/k);
  // this will throw an error if k > 2*nsamp; fix later
  for(int i = 0; i < k; i++){
    if (i < nsamp){
      mu[i] = y[step*i];
  }

  }

  double sig_ub = sig_h;

  NumericVector sig = rep(0.9*bw, mu_length);
  NumericVector sig_adj(mu_length);
  NumericVector mu_adj(mu_length);
  NumericVector p(mu_length);

  NumericVector siginit = rep(1.5*bw, nsamp);

  NumericVector loglike_curr_init(nsamp);
  for(int i = 0; i < nsamp; i++){
    loglike_curr_init[i] = KDE_density_rcpp(head(mu, k), head(sig, k), head(w,k), y[i]);
  }


  double prior_curr;
  double prior_prop;



  double adj;
  double sig_new;
  NumericVector like_prop(nsamp);
  NumericVector mu_new;



  for(int iter = 0; iter < mc; iter++){

    NumericVector drop_add = sample(NumericVector::create(-1, 0, 1), 1, false, drop_add_prob);
    double k_prop = abs(k + drop_add[0]);
    IntegerVector ind_prop = sample(k, 1, false);

    adj = 0;

    if (k_prop > k){

      mu_new = rnormmix_rcpp(y, siginit[0]);
    //  sig_new = R::rgamma(sig_shape,1/sig_rate);

   //  Rcout << sig_new << " ";

      sig_new = R::rbeta(sig_a,sig_b)*sig_ub;

    //  Rcout << sig_new << " ";

      // sig_new = std::min(sig_new,0.8*sig_h);

      adj  = -KDE_density_rcpp(y, siginit, head(w,nsamp),mu_new[0]) - R::dbeta(sig_new/sig_ub, sig_a, sig_b, true);

      for(int i = 0; i < nsamp; i++){
        like_prop[i] = R::dnorm(y[i], mu_new[0], sig_new, false);
      }

    }

    if (k_prop < k){

      adj  = KDE_density_rcpp(y, siginit, head(w,nsamp), mu[ind_prop[0]-1]) + R::dbeta(sig_new/sig_ub, sig_a, sig_b, true);

      for(int i = 0; i < nsamp; i++){
        like_prop[i] = - R::dnorm(y[i], mu[ind_prop[0]-1], sig[ind_prop[0]-1], false);
      }

    }

    if (k_prop == k){



      mu_new[0] = mu[ind_prop[0]-1] + mu_step_size*sig[ind_prop[0]-1]*(R::rnorm(0, 1));


   sig_new = sig[ind_prop[0]-1] + 0.01*(R::rnorm(0, 1));

      if (sig_new > 0.0){

        sig_new = fmod(sig_new,sig_ub);

      }

      if (sig_new < 0.0){

        sig_new = sig_ub + fmod(sig_new,sig_ub);


      }


    //  sig_new = exp(log(sig[ind_prop[0]-1]) + 0.01*(R::rnorm(0, 1)));

    // sig_new = std::min(sig_new,0.99*sig_h);






      adj = 0;

      for(int i = 0; i < nsamp; i++){
        like_prop[i] = R::dnorm(y[i], mu_new[0], sig_new, false) - R::dnorm(y[i], mu[ind_prop[0]-1], sig[ind_prop[0]-1], false);
      }

    }

    double loglike_curr = sum(loglike_curr_init);
    NumericVector loglike_prop_a = log((k*(exp(loglike_curr_init))+ like_prop)/k_prop);
    double loglike_prop = sum(loglike_prop_a);

    prior_curr = R::dgamma(k, 3.0*2.0*log(nsamp), 1.0/2.0, true);
    prior_prop = R::dgamma(k_prop, 3.0*2.0*log(nsamp), 1.0/2.0, true);

    double ratio = exp(loglike_prop - loglike_curr + prior_prop - prior_curr + adj);


    if(R::runif(0,1) < ratio){


      if (k_prop > k){

        mu[k] = mu_new[0];
        sig[k] = sig_new;
      }

      if (k == k_prop){
        mu[ind_prop[0]-1] = mu_new[0];
        sig[ind_prop[0]-1] = sig_new;
      }

      if (k_prop < k){

        mu[ind_prop[0]-1] = mu[k-1];
        sig[ind_prop[0]-1] = sig[k-1];

      }


      loglike_curr_init = loglike_prop_a;
      k = k_prop;

    }


    ksamps[iter] = k;

    if (mu_h != -100){
    for(int i = 0; i < k; i++){

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  sig_adj[i] = 1.0/pow((1.0/pow(sig[i],2) - 1.0/pow(sig_h,2)),0.5);
    mu_adj[i] = (mu[i]/pow(sig[i],2) - mu_h/pow(sig_h,2))*pow(sig_adj[i],2);
    double a = R::dnorm(mu_h, mu_adj[i], pow(pow(sig_adj[i],2) + pow(sig_h,2),0.5),false);
    p[i] = 1.0/a;

    }


   double normalizer = sum(head(p,k));
    for (int i = 0; i < k; i++){
      p[i] = k*p[i]/normalizer;

    }

    for(int i = 0; i < ygrid.length(); i++){
      fsamps(i,iter) = exp(KDE_density_rcpp(head(mu, k), head(sig, k), head(w,k), ygrid[i]));

      fsamps_adj(i,iter) = exp(KDE_density_rcpp(head(mu_adj, k), head(sig_adj, k), head(p,k), ygrid[i]));

    }

    }else{


    for(int i = 0; i < ygrid.length(); i++){
      fsamps(i,iter) = exp(KDE_density_rcpp(head(mu, k), head(sig, k), head(w,k), ygrid[i]));
      fsamps_adj(i,iter) = 1.0;
      //fsamps_adj(i,iter) = exp(KDE_density_rcpp(head(mu_adj, k), head(sig_adj, k), head(p,k), ygrid[i]));

    }

    }


  }

  return(List::create(Named("fsamps") = fsamps , Named("fsamps_adj") = fsamps_adj, Named("ksamps") = ksamps, Named("ygrid") = ygrid));
}
