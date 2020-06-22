set.seed(123)
x <- list(n = 25, r = 15, alpha = 0.0439597032152844, beta = 0.189436023504946,
          p0 = 0.45, pa = 0.7, alpha_param = 0.05, beta_param = 0.2)


expect_equivalent(
  current = as.list(fleming1stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2)),
  target = x)
