set.seed(123)
x <- list(n = 90, n_E = 45, n_C = 45, r = 7, d = 17.8, alpha = 0.0548844698916037,
     beta = 0.204267487005286, p0 = 0.45, pa = 0.7, alpha_param = 0.05,
     beta_param = 0.2)
x <- data.frame(x)
expect_equivalent(
  current = exact1stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2,alloc=1,type="normal"),
  target = x)
