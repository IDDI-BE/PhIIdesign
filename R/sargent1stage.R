
#' @title The Sargent 1-stage function
#' @description The goal of a phase II trial is to make a preliminary determination regarding the activity and
#' tolerability of a new treatment and thus to determine whether the treatment warrants
#' further study in the phase III setting. \cr
#' This function calculates the sample size needed in a Sargent 1-stage design which is a
#' three-outcome design that allows for three outcomes: reject \eqn{H(0)}, reject \eqn{H(a)}, or reject neither.
#' @param p0 probability of the uninteresting response (null hypothesis \eqn{H0})
#' @param pa probability of the interesting response (alternative hypothesis Ha)
#' @param alpha Type I error rate \eqn{P(reject H0|H0)}
#' @param beta Type II error rate \eqn{P(reject Ha|Ha)}
#' @param eta \eqn{P(reject Ha|H0)}
#' @param pi \eqn{P(reject H0|Ha)}
#' @param eps tolerance default value = 0.005
#' @param N_min minimum sample size value for grid search
#' @param N_max maximum sample size value for grid search
#' @param ... further arguments passed on to the methods
#' @return a data.frame with elements
#' \itemize{
#' \item N: total number of patients
#' \item r: cutoff point r. Note if \code{n <= r} --> futility.
#' \item s: cutoff point s. Note if \code{n >= s} --> efficacy.
#' \item eff: r/N
#' \item 90%CI_low: exact 1-2*alpha confidence interval lower limit
#' \item 90%CI_high: exact 1-2*alpha confidence interval upper limit
#' \item alpha: the actual alpha value which is smaller than \code{alpha_param + eps}
#' \item beta: the actual beta value where which is smaller than \code{beta_param + eps}
#' \item eta: the actual eta value which is smaller than \code{eta_param - eps}
#' \item pi: the actual pi value which is smaller than \code{pi_param - eps}
#' \item p0: your provided \code{p0} value
#' \item pa: your provided \code{pa} value
#' \item alpha_param: your provided \code{alpha} value
#' \item beta_param: your provided \code{beta} value
#' \item eta_param: your provided \code{eta} value
#' \item pi_param: your provided \code{pi} value
#' }
#' @references Sargent DJ, Chan V, Goldberg RM. A three-outcome design for phase II clinical trials. Control Clin Trials. 2001;22(2):117-125. doi:10.1016/s0197-2456(00)00115-x
#' @export
#' @examples
#' sargent1stage(p0 = 0.5, pa = 0.65, alpha = 0.1, beta = 0.1, eta = 0.8, pi = 0.8,
#'               eps = 0.005, N_min = 0, N_max = 100)
#' sargent1stage(p0 = 0.2, pa = 0.35, alpha = 0.1, beta = 0.1, eta = 0.8, pi = 0.8,
#'               eps = 0.005, N_min = 35, N_max = 50)
#'
#' test <- data.frame(p0 = c(0.05,0.1,0.2,0.3,0.4,0.5),
#'                    pa = c(0.05,0.1,0.2,0.3,0.4,0.5) + 0.15)
#' test <- merge(test,
#'               expand.grid(alpha = c(0.1,0.05), beta = 0.1, eta = 0.8, pi = 0.8))
#' samplesize <- sargent1stage(p0 = test$p0, pa = test$pa,
#'                             alpha = test$alpha, beta = test$beta,
#'                             eta = test$eta, pi = test$pi,
#'                             eps = 0.005, N_min = 20, N_max = 70)
#' samplesize <- lapply(seq_len(nrow(test)), FUN=function(i){
#'   setting <- test[i, ]
#'   sargent1stage(p0 = setting$p0, pa = setting$pa,
#'                 alpha = setting$alpha, beta = setting$beta, eta = setting$eta, pi = setting$pi,
#'                 eps = 0.005, N_min = 20, N_max = 70)
#' })
#' samplesize <- do.call(rbind, samplesize)
sargent1stage <- function(p0, pa, alpha, beta, eta, pi, eps = 0.005, N_min, N_max, ...) {
  if(length(p0) > 1 && length(pa) > 1){
    ## Use Rcpp implementation
    results <- mapply(null = p0, alternative = pa, alpha = alpha, beta = beta, eta = eta, pi = pi,
                      FUN = function(null, alternative, alpha, beta, eta, pi, eps, N_min, N_max){
                        sargent1stage.default(p0 = null, pa = alternative, alpha = alpha, beta = beta, eta = eta, pi = pi, eps = eps, N_min = N_min, N_max = N_max, rcpp = TRUE)
                      }, eps = eps, N_min = N_min, N_max = N_max,
                      SIMPLIFY = FALSE)
    results <- data.table::rbindlist(results)
    results <- data.table::setDF(results)
  }else{
    ## Use plain R implementation
    results <- sargent1stage.default(p0 = p0, pa = pa, alpha = alpha, beta = beta, eta = eta, pi = pi, eps = eps, N_min = N_min, N_max = N_max, ...)
  }
  results
}

sargent1stage.default <- function(p0, pa, alpha, beta, eta, pi, eps = 0.005, N_min, N_max, rcpp = TRUE) {
  # R CMD check happiness
  N <- NULL
  if (pa < p0) {
    stop("p0 should be smaller than pa")
  }

  # Get for all N's, all possible r's (0:N-2) and select where beta_temp=P(X<=r|Ha)<=beta+eps
  #------------------------------------------------------------------------------------------
  if(rcpp){
    res0 <- sargent1stage_N_r(N_min, N_max)
    res0 <- cbind(N = res0$N, r = res0$r)
    res0 <- data.frame(res0)
  }else{
    res0 <- lapply(N_min:N_max, FUN = function(a) cbind(N = a, r = 0:(a - 2)))
    res0 <- do.call(rbind, res0)
    res0 <- data.frame(res0)
  }

  res0$beta_temp <- pbinom(q = res0$r, size = res0$N, prob = pa, lower.tail = T)
  res1 <- res0[res0$beta_temp <= (beta + eps), ]
  names(res1) <- c("N", "r", "beta_temp")

  # Calculate eta for all scenarios (only r needed)
  #------------------------------------------------
  res1$eta_temp <- pbinom(q = res1$r, size = res1$N, prob = p0, lower.tail = T)
  res1$diff_eta <- eta - res1$eta_temp
  res2 <- res1[res1$diff_eta <= eps, ]

  # Get for selected N's, r's: s's ((r+2):n), and calculate pi
  #----------------------------------------------------------
  if(rcpp){
    res3 <- sargent1stage_N_r_s(N = res2$N, r = res2$r, beta_temp = res2$beta_temp, eta_temp = res2$eta_temp)
    res3 <- cbind(N = res3$N, r = res3$r, beta_temp = res3$beta_temp, eta_temp = res3$eta_temp, s = res3$s)
    res3 <- data.frame(res3)
  }else{
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

  # Calculate pi all scenarios (s needed)
  #--------------------------------------
  res3$pi_temp <- 1 - pbinom(q = res3$s - 1, size = res3$N, prob = pa, lower.tail = T)
  res3$diff_pi <- pi - res3$pi_temp
  res4 <- res3[res3$diff_pi <= eps, ]

  # Calculate alpha for all scenarios (s needed)
  #--------------------------------------------
  res4$alpha_temp <- 1 - pbinom(q = res4$s - 1, size = res4$N, prob = p0, lower.tail = T)
  res4$diff_alpha <- res4$alpha_temp - alpha
  res5 <- res4[res4$diff_alpha <= eps, ]

  res5$alpha <- alpha
  res5$beta <- beta
  res5$eta <- eta
  res5$pi <- pi

  res5 <- data.table::setDT(res5)
  res5 <- res5[, N_min := min(N)]
  res <- res5[which(N == N_min), ]
  res <- setDF(res)
  res <- data.table::setnames(res,
                              old = c("alpha", "beta", "eta", "pi"),
                              new = c("alpha_param", "beta_param", "eta_param", "pi_param"))
  res <- data.table::setnames(res,
                              old = c("alpha_temp", "beta_temp", "eta_temp", "pi_temp"),
                              new = c("alpha", "beta", "eta", "pi"))
  res <- cbind(design_nr=1:dim(res)[1],
               res[, c("r", "s", "N", "alpha", "beta", "eta", "pi")],
               p0 = p0,
               pa = pa,
               res[, c("alpha_param", "beta_param", "eta_param", "pi_param")])

  # Calculate exact 1-2*alpha confidence interval

  res$eff <- paste0(res$s, "/", res$N, " (", 100 * round((res$s) / res$N, 3), "%)")
  CI <- mapply(a = res$s, b = res$N, FUN = function(a, b) binom::binom.confint(a,b,conf.level=1-2*alpha,methods="exact"))
  res$CI_low  <- round(100 * unlist(CI[rownames(CI) == "lower", ]),2)
  res$CI_high <- round(100 * unlist(CI[rownames(CI) == "upper", ]),2)

  res<-res[, c("design_nr","N","r","s","eff","CI_low","CI_high","alpha","beta","eta","pi","p0","pa","alpha_param","beta_param","eta_param","pi_param")]

  class(res) <- c("sargent", "data.frame")
  res
}

# TEST (Sargent DJ, Goldberg RM. A Three-Outcome Design for Phase II Clinical Trials. Controlled Clinical Trials 22:117-125
#--------------------------------------------------------------------------------------------------------------------------
# test_0<- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,2),pa=rep((a+0.15),2),alpha=c(0.1,0.05),beta=c(0.1,0.1),eta=c(0.8,0.8),pi=c(0.8,0.8)),a=c(0.05,0.1,0.2,0.3,0.4,0.5),SIMPLIFY=F)))
#
# for (i in 1:dim(test_0)[1]){
#  res<- PhIIdesign::sargent1stage(p0=test_0[i,]$p0,pa=test_0[i,]$pa,alpha=test_0[i,]$alpha,beta=test_0[i,]$beta,eta=test_0[i,]$eta,pi=test_0[i,]$pi,eps=0.005,N_min=20,N_max=70)
#
#  if (i==1) {test_list      <-list(res)}   # Create list with first dataset
#  if (i!=1) {test_list[[i]] <-res }        # Next iterations: append dataset
#
# }
# test<- data.frame(do.call("rbind",test_list))
# test<-test[,c("p0","pa","alpha_param","beta_param","eta_param","pi_param","alpha","beta","eta","pi","r","s","N")]
