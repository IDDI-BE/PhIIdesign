set.seed(123)
x <- list(r = 32, s = 34, N = 57, alpha = 0.0924252412082096, beta = 0.104393155688326,
          eta = 0.855378126445489, pi = 0.838046143845541, p0 = 0.5,
          pa = 0.65, alpha_param = 0.1, beta_param = 0.1, eta_param = 0.8,
          pi_param = 0.8)
expect_equivalent(
  current = as.list(sargent1stage(p0=0.5,pa=0.65,alpha=0.1,beta=0.1,eta=0.8,pi=0.8,eps = 0.005,N_min=0,N_max=100)),
  target = x)
