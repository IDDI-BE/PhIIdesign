set.seed(123)
x <- list(design_nr=c(1:3),
          r1 = c(1, 1, 1),
          n1 = c(14, 14, 12),
          r2 = c(3, 3, 3),
          s = c(5, 6, 6),
          n2 = c(7, 11, 14),
          N = c(21, 25, 26),
          eff = c("5/21 (23.8%)", "6/25 (24%)", "6/26 (23.1%)"),
          "90%CI_low" = c(9.94, 11.09, 10.83),
          "90%CI_high" = c(38.71, 38, 38.15),
          EN.p0 = c(16.9075960163903, 18.5690794543276, 16.773968474954),
          PET.p0 = c(0.58462914051567, 0.58462914051567, 0.659002251789),
          MIN = c("Minimax", "", ""),
          OPT = c("", "", "Optimal"),
          ADMISS = c("", "", ""),
          alpha = c(0.0511412111636121, 0.0323533353291087, 0.0359671496394091),
          beta = c(0.100822444566976, 0.0641267609811352, 0.0946167011418047),
          eta = c(0.857801376756208, 0.799703105552571, 0.81304011630856),
          pi = c(0.796294384317041, 0.797423012704184, 0.804780447494023),
          lambda = c(0.0910574120801796, 0.167943559118321, 0.150992734052031),
          delta = c(0.102883171115982, 0.138450226314681, 0.100602851364173),
          p0 = c(0.1, 0.1, 0.1),
          pa = c(0.3, 0.3, 0.3),
          alpha_param = c(0.05, 0.05, 0.05),
          beta_param = c(0.1, 0.1, 0.1),
          eta_param = c(0.8, 0.8, 0.8),
          pi_param = c(0.8, 0.8, 0.8))

expect_equivalent(
  current = as.list(sargent2stage(p0=0.1,pa=0.3,alpha=0.05,beta=0.1,eta=0.8,pi=0.8,eps = 0.005,N_min=15,N_max=30)),
  target = x)

expect_equivalent(
  current = as.list(sargent2stage(p0=0.1,pa=0.3,alpha=0.05,beta=0.1,eta=0.8,pi=0.8,eps = 0.005,N_min=15,N_max=30, method = "original")),
  target = as.list(sargent2stage(p0=0.1,pa=0.3,alpha=0.05,beta=0.1,eta=0.8,pi=0.8,eps = 0.005,N_min=15,N_max=30, method = "speedup")))

expect_equal(
  current = PhIIdesign:::probha(n1 = 14, n2 = 7, r1 = 1, s = 5, p = 0.1, type = "original"),
  target = PhIIdesign:::probha(n1 = 14, n2 = 7, r1 = 1, s = 5, p = 0.1, type = "speedup"))
