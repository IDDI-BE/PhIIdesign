set.seed(123)
x <- list(design_nr=1, N = 25, r = 15, alpha = 0.0439597032152844, beta = 0.189436023504946,
          p0 = 0.45, pa = 0.7, alpha_param = 0.05, beta_param = 0.2)


expect_equivalent(
  current = as.list(fleming1stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2))[c("design_nr","N","r","alpha","beta","p0","pa","alpha_param","beta_param")],
  target = x)


expect_equal(
  as.list(fleming1stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2)[,c("N","r","alpha","beta","p0","pa","alpha_param","beta_param")]),
  PhIIdesign:::fleming_single_stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2)
)
