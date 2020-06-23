#include <Rcpp.h>


// [[Rcpp::export]]
Rcpp::List fleming_single_stage(double p0, double pa, double alpha = 0.05, double beta = 0.2, double eps = 0.005) {
  int n = 0;
  double beta_temp = 1;
  double r_temp = -1;
  while (beta_temp > beta + eps) {
    if(n % 1000 == 0){
      Rcpp::checkUserInterrupt();
    }
    n = n + 1;
    r_temp = R::qbinom(1 - (alpha + eps), n, p0, true, false);
    beta_temp = R::pbinom(r_temp, n, pa, true, false);
  }
  double alpha_temp = 1 - R::pbinom(r_temp, n, p0, true, false);

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("n") = n,
    Rcpp::Named("r") = r_temp,
    Rcpp::Named("alpha") = alpha_temp,
    Rcpp::Named("beta") = beta_temp,
    Rcpp::Named("p0") = p0,
    Rcpp::Named("pa") = pa,
    Rcpp::Named("alpha_param") = alpha,
    Rcpp::Named("beta_param") = beta);
  return output;
}


