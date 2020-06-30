
#----------------------------------------------------------------------------------------#
# Cumulative distribution function for difference between two binomial variables Z=X1-X2 #
# (allowing for unequal sample sizes in both groups)                                     #
#----------------------------------------------------------------------------------------#

bin_dif_cdf <- function(z, n1, n2, p1, p2, type = c("exact", "normal")) { # z=cutoff difference between two binomial variables; Z=X1-X2
  type <- match.arg(type)
  if (type == "exact") {
    prob <- 0

    for (diff in -n2:z) { # Z=X1-X2, so -n2 is minimum
      i <- (0:max(n1, n2))
      prob <- prob + sum(dbinom(x = i + diff, size = n1, prob = p1) * dbinom(x = i, size = n2, prob = p2))
    }
  }else if (type == "normal") {
    prob <- pnorm(q = z + 0.5, mean = (n1 * p1 - n2 * p2), sd = (n1 * p1 * (1 - p1) + n2 * p2 * (1 - p2))**0.5, lower.tail = T)
  }
  prob
}

#' @title The exact 1-stage function, randomized
#' @description This function calculates sample sizes of an exact randomized 1-stage design, for p0<pa
#' @param p0 probability of the uninteresting response (null hypothesis H0)
#' @param pa probability of the interesting response (alternative hypothesis Ha)
#' @param alpha Type I error rate \eqn{P(reject H0|H0)}
#' @param beta Type II error rate \eqn{P(reject Ha|Ha)}
#' @param eps tolerance default value = 0.005
#' @param alloc allocation ratio (e.g. 2 for 2:1). alloc is assumed >=1 (more patients in the experimental arm)
#' @param ... refers to type="exact" or "normal" for pdf and cdf difference of two binomial variables
#' @return a data.frame with elements
#' \itemize{
#' \item n: total number of patients
#' \item n_E: number of patients in experimental arm
#' \item n_C: number of patients in control arm
#' \item r: TODO. Variable of interest: diff=E-C. Note if (E-C)<=r --> futility\cr
#' \item d: TODO. Variable of interest: diff=E-C. Note if if (E-C)> r --> efficacy
#' \item alpha: the actual alpha value which is smaller than \code{alpha_param + eps}
#' \item beta: the actual beta value where which is smaller than \code{beta_param + eps}
#' \item p0: your provided \code{p0} value
#' \item pa: your provided \code{pa} value
#' \item alpha_param: your provided \code{alpha} value
#' \item beta_param: your provided \code{beta} value
#' }
#' @details
#' Note that all sample sizes can be entered, but only calculations for sample sizes, without any decimals (taking into account allocation ratio)
#' @export
#' @examples
#' exact1stage(p0 = 0.45, pa = 0.7, alpha = 0.05, beta = 0.2, alloc = 1, type = "normal")
#' exact1stage(p0 = 0.45, pa = 0.7, alpha = 0.05, beta = 0.2, alloc = 1, type = "exact")
exact1stage <- function(p0, pa, alpha, beta, eps = 0.005, alloc = 1, ...) {
  stopifnot(alloc >= 1)
  if (pa < p0) {
    stop("p0 should be smaller than pa")
  }

  n_C <- 0
  beta_temp <- 1

  while (beta_temp > beta + eps) {
    n_C <- n_C + 1
    n_E <- n_C * alloc

    r_temp <- -1
    alpha_temp <- 1

    while (alpha_temp > alpha + eps & r_temp <= n_E - 1) {
      r_temp <- r_temp + 1
      alpha_temp <- 1 - bin_dif_cdf(z = r_temp, n1 = n_E, n2 = n_C, p1 = (p0 + pa) / 2, p2 = (p0 + pa) / 2, ...)
    }
    beta_temp <- bin_dif_cdf(z = r_temp, n1 = n_E, n2 = n_C, p1 = pa, p2 = p0, ...)
  }

  n <- n_E + n_C

  if (alloc == 1) {
    d <- round(100 * (r_temp + 1) / n_E, 1)
    res <- data.frame(n = n, n_E = n_E, n_C = n_C, r = r_temp, d = d, alpha = alpha_temp, beta = beta_temp, p0 = p0, pa = pa, alpha_param = alpha, beta_param = beta)
  }else if (alloc > 1) {
    res <- data.frame(n = n, n_E = n_E, n_C = n_C, r = r_temp, d = NA_real_, alpha = alpha_temp, beta = beta_temp, p0 = p0, pa = pa, alpha_param = alpha, beta_param = beta)
  }
  return(res)
}

# TEST (Sargent DJ, Goldberg RM. A Three-Outcome Design for Phase II Clinical Trials. Controlled Clinical Trials 22:117-125). Table I (Two-outcome designs)
#----------------------------------------------------------------------------------------------------------------------------------------------------------
# test_0_a<- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,2),pa=rep((a+0.2 ),2),alpha=c(0.05,0.1),beta=c(0.1,0.1)),a=c(0.1,0.2,0.3,0.4),SIMPLIFY=F)))
# test_0_b<- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,2),pa=rep((a+0.15),2),alpha=c(0.05,0.1),beta=c(0.1,0.1)),a=c(0.1,0.2,0.3,0.4),SIMPLIFY=F)))
# test_0<-rbind(test_0_a,test_0_b)
#
# # Normal approximation (with continuity correction)
# for (i in 1:dim(test_0)[1]){
#
#   res<- PhIIdesign::exact1stage(p0=test_0[i,]$p0,pa=test_0[i,]$pa,alpha=test_0[i,]$alpha,beta=test_0[i,]$beta,eps=0.005,alloc=1,type="normal")
#   if (i==1) {test_list      <-list(res)}   # Create list with first dataset
#   if (i!=1) {test_list[[i]] <-res      }        # Next iterations: append dataset
#
# }
# test_normal<- data.frame(do.call("rbind",test_list))
#
# # Exact binomial test
# for (i in 1:dim(test_0)[1]){
#
#   res<- PhIIdesign::exact1stage(p0=test_0[i,]$p0,pa=test_0[i,]$pa,alpha=test_0[i,]$alpha,beta=test_0[i,]$beta,eps=0.005,alloc=1,type="exact")
#   if (i==1) {test_list      <-list(res)}   # Create list with first dataset
#   if (i!=1) {test_list[[i]] <-res      }        # Next iterations: append dataset
#
# }
# test_exact<- data.frame(do.call("rbind",test_list))
