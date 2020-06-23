
#' @title The Fleming 1-stage function
#' @description Therapeutic efficacy in clinical trials is often evaluated principally on the basis of the probability p that
#' an eligible patient receiving the treatment regimen will experience a regression (like a tumor e.g.).\cr
#' This function calculates sample sizes of the Fleming single-stage design, for \eqn{p_0 < p_a} for the requested Type I (alpha) and Type II error (beta).
#' @param p0 probability of the uninteresting response (null hypothesis H0)
#' @param pa probability of the interesting response (alternative hypothesis Ha)
#' @param alpha Type I error rate \eqn{P(reject H0|H0)}
#' @param beta Type II error rate \eqn{P(reject Ha|Ha)}
#' @param eps tolerance  default value = 0.005
#' @details
#' TODO:
#' actual eta>=eta-eps\cr
#' actual p>=pi-eps
#' @return a data.frame with elements
#' \itemize{
#' \item n: total number of patients
#' \item r: quantile function of \code{1 - (alpha + eps)} at \code{n} under \code{p0}. Note if \code{n <= r} --> futility.
#' \item alpha: the actual alpha value where your provided \code{alpha_param <= alpha + eps}
#' \item beta: the actual beta value where your provided \code{beta_param <= beta + eps}
#' \item p0: your provided \code{p0} value
#' \item pa: your provided \code{pa} value
#' \item alpha_param: your provided \code{alpha} value
#' \item beta_param: your provided \code{beta} value
#' }
#' @references Fleming TR. One-sample multiple testing procedure for phase II clinical trials. Biometrics. 1982;38(1):143-151.
#' @export
#' @examples
#' fleming1stage(p0 = 0.45, pa = 0.7, alpha = 0.05, beta = 0.2)
#' fleming1stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.1, eps = 0.005)
#' fleming1stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.1, eps = 0.00001)
#'
#' ## For several combinations of p0 and pa
#' ## not it is important that p0 is not equal to pa
#' test <- expand.grid(p0 = seq(0, 0.95, by = 0.05),
#'                     pa = seq(0, 0.95, by = 0.05))
#' test <- subset(test, (pa - p0) > 0.00001)
#' samplesize <- fleming1stage(p0 = test$p0, pa = test$pa, alpha = 0.05, beta = 0.2, eps = 0.0005)
#' samplesize <- fleming1stage(p0 = test$p0, pa = test$pa, alpha = 0.05, beta = 0.1, eps = 0.0005)
#' samplesize <- fleming1stage(p0 = test$p0, pa = test$pa, alpha = 0.01, beta = 0.2, eps = 0.0005)
#' samplesize <- fleming1stage(p0 = test$p0, pa = test$pa, alpha = 0.01, beta = 0.1, eps = 0.0005)
#'
#' ## these 2 are the same
#' samplesize <- fleming1stage(p0 = test$p0, pa = test$pa, alpha = 0.05, beta = 0.2)
#' samplesize <- mapply(p0 = test$p0, pa = test$pa, FUN=function(p0, pa){
#'   fleming1stage(p0 = p0, pa = pa, alpha = 0.05, beta = 0.2)
#' }, SIMPLIFY = FALSE)
#' samplesize <- do.call(rbind, samplesize)
fleming1stage <- function(p0, pa, alpha = 0.05, beta = 0.2, eps = 0.005){
  stopifnot(length(alpha) == 1 && length(beta) == 1 && length(eps) == 1)
  stopifnot(all(p0 >= 0) && all(p0 <= 1))
  stopifnot(all(pa >= 0) && all(pa <= 1))
  stopifnot(alpha >= 0 && alpha <= 1)
  stopifnot(beta >= 0 && beta <= 1)
  stopifnot((beta + eps) <= 1)
  if (!all(p0 < pa)) {
    stop("p0 should be smaller than pa")
  }
  if(length(p0) > 1 && length(pa) > 1){
    ## Use Rcpp implementation
    results <- mapply(null = p0, alternative = pa,
                      FUN = function(null, alternative, alpha, beta, eps){
                        fleming_single_stage(p0 = null, pa = alternative, alpha = alpha, beta = beta, eps = eps)
                        }, alpha = alpha, beta = beta, eps = eps,
                      SIMPLIFY = FALSE)
    results <- data.table::rbindlist(results)
    results <- data.table::setDF(results)
  }else{
    ## Use plain R implementation
    results <- fleming1stage.default(p0 = p0, pa = pa, alpha = alpha, beta = beta, eps = eps)
  }
  results
}

fleming1stage.default <- function(p0, pa, alpha = 0.05, beta = 0.2, eps = 0.005) {
  stopifnot(p0 >= 0 && p0 <= 1)
  stopifnot(pa >= 0 && pa <= 1)
  stopifnot(alpha >= 0 && alpha <= 1)
  stopifnot(beta >= 0 && beta <= 1)
  stopifnot((beta + eps) <= 1)
  if(!(p0 < pa)) {
    stop("p0 should be smaller than pa")
  }

  n <- 0
  beta_temp <- 1

  while (beta_temp > beta + eps) {
    n <- n + 1
    r_temp <- qbinom(p = 1 - (alpha + eps), size = n, prob = p0, lower.tail = T)
    beta_temp <- pbinom(q = r_temp, size = n, prob = pa, lower.tail = T)
  }

  alpha_temp <- 1 - pbinom(q = r_temp, size = n, prob = p0, lower.tail = T)
  res <- data.frame(n = n,
                    r = r_temp,
                    alpha = alpha_temp,
                    beta = beta_temp,
                    p0 = p0,
                    pa = pa,
                    alpha_param = alpha,
                    beta_param = beta)
  res
}

# TEST (A'Hern RP. Sample size tables for exact single-stage phase II designs. Statistics in Medicine 2001;20:859-866: Table 1
#------------------------------------------------------------------------------------------------------------------------------
# test0                    <- data.frame(do.call("rbind", mapply(function(a) cbind(p0=a,pa=seq(a,0.95,by=0.05)),a=seq(0.05,0.95,by=0.05),SIMPLIFY=F)))
# test                     <- test0[test0$pa>test0$p0,]
# test_alpha_0.05_beta_0.2 <- cbind(test,data.frame(do.call("rbind",mapply(function (a,b) PhIIdesign::fleming1stage(p0=a,pa=b,alpha=0.05,beta=0.2,eps=0.0005),a=test$p0,b=test$pa,SIMPLIFY=F))))
# test_alpha_0.05_beta_0.1 <- cbind(test,data.frame(do.call("rbind",mapply(function (a,b) PhIIdesign::fleming1stage(p0=a,pa=b,alpha=0.05,beta=0.1,eps=0.0005),a=test$p0,b=test$pa,SIMPLIFY=F))))
# test_alpha_0.01_beta_0.2 <- cbind(test,data.frame(do.call("rbind",mapply(function (a,b) PhIIdesign::fleming1stage(p0=a,pa=b,alpha=0.01,beta=0.2,eps=0.0005),a=test$p0,b=test$pa,SIMPLIFY=F))))
# test_alpha_0.01_beta_0.1 <- cbind(test,data.frame(do.call("rbind",mapply(function (a,b) PhIIdesign::fleming1stage(p0=a,pa=b,alpha=0.01,beta=0.1,eps=0.0005),a=test$p0,b=test$pa,SIMPLIFY=F))))
