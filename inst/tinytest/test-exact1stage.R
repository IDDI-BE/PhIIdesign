set.seed(123)
x <- list(n = 90, n_E = 45, n_C = 45, r = 7, d = 17.8, alpha = 0.0548844698916037,
     beta = 0.204267487005286, p0 = 0.45, pa = 0.7, alpha_param = 0.05,
     beta_param = 0.2)
x <- data.frame(x)
expect_equivalent(
  current = exact1stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2,alloc=1,type="normal"),
  target = x)
expect_equivalent(
  current = exact1stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2,alloc=1,type="normal"),
  target = data.frame(PhIIdesign:::rcpp_exact1stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2,alloc=1,type="normal")))
expect_equivalent(
  current = exact1stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2,alloc=1,type="exact"),
  target = data.frame(PhIIdesign:::rcpp_exact1stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2,alloc=1,type="exact")))

test <- rbind(
  data.frame(p0 = c(0.1, 0.2, 0.3, 0.4),
             pa = c(0.1, 0.2, 0.3, 0.4) + 0.2),
  data.frame(p0 = c(0.1, 0.2, 0.3, 0.4),
             pa = c(0.1, 0.2, 0.3, 0.4) + 0.15))
test <- merge(test,
              expand.grid(alpha = c(0.05, 0.1), beta = 0.1))

expect_equivalent(
  current = do.call(rbind, mapply(p0 = test$p0, pa = test$pa, alpha = test$alpha, beta = test$beta,
                                  FUN=function(p0, pa, alpha, beta){
                                    exact1stage(p0 = p0, pa = pa, alpha = alpha, beta = beta,
                                                eps = 0.005, alloc = 1, type = "normal")
                                  }, SIMPLIFY = FALSE)),
  target = exact1stage(p0 = test$p0, pa = test$pa,
                       alpha = test$alpha, beta = test$beta,
                       eps = 0.005, alloc = 1, type = "normal"))
expect_equivalent(
  current = do.call(rbind, mapply(p0 = test$p0, pa = test$pa, alpha = test$alpha, beta = test$beta,
                                  FUN=function(p0, pa, alpha, beta){
                                    exact1stage(p0 = p0, pa = pa, alpha = alpha, beta = beta,
                                                eps = 0.1, alloc = 1, type = "exact")
                                  }, SIMPLIFY = FALSE)),
  target = exact1stage(p0 = test$p0, pa = test$pa,
                       alpha = test$alpha, beta = test$beta,
                       eps = 0.1, alloc = 1, type = "exact"))


expect_equivalent(
  current = PhIIdesign:::bin_dif_cdf(0.306526, 0, 2, 0.575, 0.575, type="normal"),
  target = PhIIdesign:::rcpp_bin_dif_cdf(0.306526, 0, 2, 0.575, 0.575, type="normal"))
expect_equivalent(
  current = PhIIdesign:::bin_dif_cdf(0.306526, 1, 5, 0.1, 0.575, type="exact"),
  target = PhIIdesign:::rcpp_bin_dif_cdf(0.306526, 1, 5, 0.1, 0.575, type="exact"))
expect_equivalent(
  current = PhIIdesign:::bin_dif_cdf(0.1, 90, 100, 0.1, 0.05, type="exact"),
  target = PhIIdesign:::rcpp_bin_dif_cdf(0.1, 90, 100, 0.1, 0.05, type="exact"))


