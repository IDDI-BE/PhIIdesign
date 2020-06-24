library(benchr)
library(ggplot2)
library(Rcpp)
N_min <- 0
N_max <- 100

cppFunction(
  'Rcpp::List sargent1stage_N_r(int N_min, int N_max) {
        std::vector<int> _N;
        std::vector<int> _r;
        for(int N = N_min; N <= N_max; N++){
          for(int r = 0; r <= N-2; r++){
            _N.push_back (N);
            _r.push_back (r);
          }
        }
        Rcpp::List output = Rcpp::List::create(
          Rcpp::Named("N") = _N,
          Rcpp::Named("r") = _r);
        return output;
    }')


comparison <- benchmark(
  {
    res0 <- sargent1stage_N_r(N_min, N_max)
    res0 <- cbind(N = res0$N, r = res0$r)
    res0 <- data.frame(res0)
  },
  {
    res0 <- data.table(N = N_min:N_max)
    res0 <- res0[, list(r = 0:(.BY$N - 2)), by = list(N)]
    res0 <- setDF(res0)
  },
  {
    res0 <- lapply(N_min:N_max, FUN = function(a) cbind(N = a, r = 0:(a - 2)))
    res0 <- do.call(rbind, res0)
    res0 <- data.frame(res0)
  }
)
comparison
boxplot(comparison) + coord_flip()

cppFunction(
  'Rcpp::List sargent1stage_N_r_s(Rcpp::IntegerVector N, Rcpp::IntegerVector r, Rcpp::NumericVector beta_temp, Rcpp::NumericVector eta_temp) {
        std::vector<int> _N;
        std::vector<int> _r;
        std::vector<double> _beta_temp;
        std::vector<double> _eta_temp;
        std::vector<int> _s;
        for(int i = 0; i < N.size(); i++){
          int Ni = N[i];
          int ri = r[i];
          for(int s = ri + 2; s <= Ni; s++){
            _N.push_back(N[i]);
            _r.push_back(r[i]);
            _beta_temp.push_back(beta_temp[i]);
            _eta_temp.push_back(eta_temp[i]);
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
    }')

N_min <- 20
N_max <- 100
eps = 0.005
alpha = 0.1
beta = 0.1
eta = 0.8
pi = 0.8
p0 = 0.2
pa = 0.35
res0 <- lapply(N_min:N_max, FUN = function(a) cbind(N = a, r = 0:(a - 2)))
res0 <- do.call(rbind, res0)
res0 <- data.frame(res0)
res0$beta_temp <- pbinom(q = res0$r, size = res0$N, prob = pa, lower.tail = T)
res1 <- res0[res0$beta_temp <= (beta + eps), ]
names(res1) <- c("N", "r", "beta_temp")

# Calculate eta for all scenarios (only r needed)
#------------------------------------------------
res1$eta_temp <- pbinom(q = res1$r, size = res1$N, prob = p0, lower.tail = T)
res1$diff_eta <- eta - res1$eta_temp
res2 <- res1[res1$diff_eta <= eps, ]



comparison <- benchmark(
  {
    res3 <- sargent1stage_N_r_s(N = res2$N, r = res2$r, beta_temp = res2$beta_temp, eta_temp = res2$eta_temp)
    res3 <- cbind(N = res3$N, r = res3$r, beta_temp = res3$beta_temp, eta_temp = res3$eta_temp, s = res3$s)
    res3 <- data.frame(res3)
  },
  {
    res3 <- setDT(res2)
    res3 <- res3[, list(s = (.BY$r + 2):.BY$N), by = list(N, r, beta_temp, eta_temp)]
    res3 <- setDF(res3)
  },
  {
    res3 <- mapply(
      FUN = function(N, r, beta_temp, eta_temp){
        cbind(N = N, r = r, beta_temp = beta_temp, eta_temp = eta_temp, s = (r + 2):N)
      },
      N = res2$N,
      r = res2$r,
      beta_temp = res2$beta_temp,
      eta_temp = res2$eta_temp)
    res3 <- do.call(rbind, res3)
    res3 <- data.frame(res3)
  }
)
comparison
boxplot(comparison) + coord_flip()
