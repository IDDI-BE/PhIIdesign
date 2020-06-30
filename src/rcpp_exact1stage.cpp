#include <Rcpp.h>


// [[Rcpp::export]]
double rcpp_bin_dif_cdf(double z, int n1, int n2, double p1, double p2, std::string type) {
  double prob = 0;
  if (type == "exact") {
    // Z=X1-X2, so -n2 is minimum
    for(int diff = -n2; diff <= z; diff++){
      Rcpp::IntegerVector i = Rcpp::seq(0, std::max(n1, n2));
      prob = prob + Rcpp::sum(Rcpp::dbinom(i + diff, n1, p1) * Rcpp::dbinom(i, n2, p2));
    }
  }else if (type == "normal") {
    prob = R::pnorm(z + 0.5, (n1 * p1 - n2 * p2), std::sqrt((n1 * p1 * (1 - p1) + n2 * p2 * (1 - p2))), true, false);
  }
  return prob;
}


// [[Rcpp::export]]
Rcpp::List rcpp_exact1stage(double p0, double pa, double alpha, double beta, double eps = 0.005, int alloc = 1, std::string type = "exact") {
  int n_E = 0;
  int r_temp = -1;
  double alpha_temp = -1;
  int n_C = 0;
  double beta_temp = 1;

  while (beta_temp > beta + eps) {
    if(n_C % 1000 == 0){
      Rcpp::checkUserInterrupt();
    }
    n_C = n_C + 1;
    n_E = n_C * alloc;

    r_temp = -1;
    alpha_temp = 1;

    while (alpha_temp > alpha + eps && r_temp <= n_E - 1) {
      if(r_temp % 1000 == 0){
        Rcpp::checkUserInterrupt();
      }
      r_temp = r_temp + 1;
      alpha_temp = 1 - rcpp_bin_dif_cdf(r_temp,n_E, n_C, (p0 + pa) / 2, (p0 + pa) / 2, type);
    }
    beta_temp = rcpp_bin_dif_cdf(r_temp, n_E, n_C, pa, p0, type);
  }

  int n = n_E + n_C;

  Rcpp::NumericVector d;
  if (alloc == 1) {
    d = Rcpp::NumericVector::create(100 * ((double)r_temp + 1) / (double)n_E);
    d = Rcpp::round(d, 1);
  }else if (alloc > 1) {
    d = Rcpp::NumericVector::create(NA_REAL);
  }

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("n") = n,
    Rcpp::Named("n_E") = n_E,
    Rcpp::Named("n_C") = n_C,
    Rcpp::Named("r") = r_temp,
    Rcpp::Named("d") = d,
    Rcpp::Named("alpha") = alpha_temp,
    Rcpp::Named("beta") = beta_temp,
    Rcpp::Named("p0") = p0,
    Rcpp::Named("pa") = pa,
    Rcpp::Named("alpha_param") = alpha,
    Rcpp::Named("beta_param") = beta);
  return output;
}
