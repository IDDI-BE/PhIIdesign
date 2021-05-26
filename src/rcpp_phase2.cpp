#include <Rcpp.h>


// [[Rcpp::export]]
Rcpp::List fleming_single_stage(double p0, double pa, double alpha = 0.05, double beta = 0.2, double eps = 0.005) {
  int N = 0;
  double beta_temp = 1;
  double r_temp = -1;
  while (beta_temp > beta + eps) {
    if(N % 1000 == 0){
      Rcpp::checkUserInterrupt();
    }
    N = N + 1;
    r_temp = R::qbinom(1 - (alpha + eps), N, p0, true, false);
    beta_temp = R::pbinom(r_temp, N, pa, true, false);
  }
  double alpha_temp = 1 - R::pbinom(r_temp, N, p0, true, false);

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("N") = N,
    Rcpp::Named("r") = r_temp,
    Rcpp::Named("alpha") = alpha_temp,
    Rcpp::Named("beta") = beta_temp,
    Rcpp::Named("p0") = p0,
    Rcpp::Named("pa") = pa,
    Rcpp::Named("alpha_param") = alpha,
    Rcpp::Named("beta_param") = beta);
  return output;
}


// [[Rcpp::export]]
Rcpp::List sargent1stage_N_r(int N_min, int N_max) {
  std::vector<int> _N;
  std::vector<int> _r;
  for(int N = N_min; N <= N_max; N++){
    if(N % 1000 == 0){
      Rcpp::checkUserInterrupt();
    }
    for(int r = 0; r <= N-2; r++){
      _N.push_back (N);
      _r.push_back (r);
    }
  }
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("N") = _N,
    Rcpp::Named("r") = _r);
  return output;
}


// [[Rcpp::export]]
Rcpp::List sargent1stage_N_r_s(Rcpp::IntegerVector N, Rcpp::IntegerVector r, Rcpp::NumericVector beta_temp, Rcpp::NumericVector eta_temp) {
  std::vector<int> _N;
  std::vector<int> _r;
  std::vector<double> _beta_temp;
  std::vector<double> _eta_temp;
  std::vector<int> _s;
  for(int i = 0; i < N.size(); i++){
    int N_i = N(i);
    int r_i = r(i);
    double beta_temp_i = beta_temp(i);
    double eta_temp_i = eta_temp(i);
    for(int s = r_i + 2; s <= N_i; s++){
      _N.push_back(N_i);
      _r.push_back(r_i);
      _beta_temp.push_back(beta_temp_i);
      _eta_temp.push_back(eta_temp_i);
      _s.push_back(s);
    }
  }
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("N") = _N,
    Rcpp::Named("r") = _r,
    Rcpp::Named("beta_temp") = _beta_temp,
    Rcpp::Named("eta_temp") = _eta_temp,
    Rcpp::Named("s") = _s);
  return output;
}
