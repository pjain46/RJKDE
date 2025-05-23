// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rnormmix_rcpp
NumericVector rnormmix_rcpp(NumericVector mu, double sigma);
RcppExport SEXP _RJKDE_rnormmix_rcpp(SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rnormmix_rcpp(mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// KDE_density_rcpp
double KDE_density_rcpp(NumericVector mu, NumericVector sigma, NumericVector w, double x);
RcppExport SEXP _RJKDE_KDE_density_rcpp(SEXP muSEXP, SEXP sigmaSEXP, SEXP wSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(KDE_density_rcpp(mu, sigma, w, x));
    return rcpp_result_gen;
END_RCPP
}
// rj_mcmc_rcpp
List rj_mcmc_rcpp(NumericVector y, NumericVector ygrid, NumericVector drop_add_prob, double sig_a, double sig_b, double mu_step_size, int k, double bw, int mc, double mu_h, double sig_h);
RcppExport SEXP _RJKDE_rj_mcmc_rcpp(SEXP ySEXP, SEXP ygridSEXP, SEXP drop_add_probSEXP, SEXP sig_aSEXP, SEXP sig_bSEXP, SEXP mu_step_sizeSEXP, SEXP kSEXP, SEXP bwSEXP, SEXP mcSEXP, SEXP mu_hSEXP, SEXP sig_hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ygrid(ygridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type drop_add_prob(drop_add_probSEXP);
    Rcpp::traits::input_parameter< double >::type sig_a(sig_aSEXP);
    Rcpp::traits::input_parameter< double >::type sig_b(sig_bSEXP);
    Rcpp::traits::input_parameter< double >::type mu_step_size(mu_step_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type mc(mcSEXP);
    Rcpp::traits::input_parameter< double >::type mu_h(mu_hSEXP);
    Rcpp::traits::input_parameter< double >::type sig_h(sig_hSEXP);
    rcpp_result_gen = Rcpp::wrap(rj_mcmc_rcpp(y, ygrid, drop_add_prob, sig_a, sig_b, mu_step_size, k, bw, mc, mu_h, sig_h));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RJKDE_rnormmix_rcpp", (DL_FUNC) &_RJKDE_rnormmix_rcpp, 2},
    {"_RJKDE_KDE_density_rcpp", (DL_FUNC) &_RJKDE_KDE_density_rcpp, 4},
    {"_RJKDE_rj_mcmc_rcpp", (DL_FUNC) &_RJKDE_rj_mcmc_rcpp, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_RJKDE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
