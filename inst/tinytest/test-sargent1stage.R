set.seed(123)
x <- list(r = 32, s = 34, N = 57, alpha = 0.0924252412082096, beta = 0.104393155688326,
          eta = 0.855378126445489, pi = 0.838046143845541, p0 = 0.5,
          pa = 0.65, alpha_param = 0.1, beta_param = 0.1, eta_param = 0.8,
          pi_param = 0.8)
expect_equivalent(
  current = as.list(sargent1stage(p0=0.5,pa=0.65,alpha=0.1,beta=0.1,eta=0.8,pi=0.8,eps = 0.005,N_min=0,N_max=100, rcpp=FALSE)[,-c(1)]), # '-c(1) to omit 'design_nr' column
  target = x)

expect_equivalent(
  current = as.list(sargent1stage(p0=0.5,pa=0.65,alpha=0.1,beta=0.1,eta=0.8,pi=0.8,eps = 0.005,N_min=0,N_max=100, rcpp=TRUE)[,-c(1)]), # '-c(1) to omit 'design_nr' column
  target = x)


test <- data.frame(p0 = c(0.05,0.1,0.2,0.3,0.4,0.5),
                   pa = c(0.05,0.1,0.2,0.3,0.4,0.5) + 0.15)
test <- merge(test,
              expand.grid(alpha = c(0.1,0.05), beta = 0.1, eta = 0.8, pi = 0.8))

expect_equivalent(
  current = do.call(rbind,
                    lapply(seq_len(nrow(test)), FUN=function(i){
                      setting <- test[i, ]
                      sargent1stage(p0 = setting$p0, pa = setting$pa,
                                    alpha = setting$alpha, beta = setting$beta, eta = setting$eta, pi = setting$pi,
                                    eps = 0.005, N_min = 20, N_max = 70)
                    })),
  target = sargent1stage(p0 = test$p0, pa = test$pa,
                         alpha = test$alpha, beta = test$beta,
                         eta = test$eta, pi = test$pi,
                         eps = 0.005, N_min = 20, N_max = 70))
