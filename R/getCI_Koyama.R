
#' @title Koyama confidence interval.
#' @description
#' Calculates the two sided 1-2*alpha confidence interval based, by inverting the "rejection probability
#' function" (or the "p value function") at every possible value for the null hypothesis between [0,1].
#' A two-sided (1 − 2*alpha) CI interval includes then all values of binomial proportions for
#' which the p value for testing H0, lays within the interval [alpha, 1 − alpha]
#' @param k overall observed responses (must be larger than r1).
#' @param r1 critical value for the first stage.
#' @param n1 sample size for the first stage.
#' @param n overall sample size.
#' @param alpha determining the two sided 1-2*alpha confidence interval.
#' @param precision gives the precision (in decimal numbers) to which the confidence interval should be calculated (should be less than 10).
#' @references Koyama T and Chen H (2008): Proper inference from Simon's two-stage designs. Statistics in Medicine, 27(16):3145-3154.
#' corresponds with the "R" approach in Chan et al.
#' @export
#' @examples
#' #Calculate a Simon's two-stage design
#' res <- simon2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.2,eps = 0.005,
#' N_min = 1, N_max = 50)
#' CI <- mapply(a = res$r2 + 1, b = res$r1, c = res$n1, d = res$N, FUN = function(a, b, c, d)
#' getCI_Koyama(k = a, r1 = b, n1 = c, n = d, alpha = 0.05, precision = 3))

getCI_Koyama<-function (k,r1,n1,n,alpha,precision=3){

  stopifnot(class(k)         == "numeric" | class(k) == "integer", k >= 0)
  stopifnot(class(r1)        == "numeric", r1 >= 0)
  stopifnot(class(n1)        == "numeric", n1 >= 0, n1 > r1)
  stopifnot(class(n)         == "numeric", n >= 0, n > n1)
  stopifnot(class(alpha)     == "numeric", alpha > 0, alpha <= 1)
  stopifnot(class(precision) == "numeric", precision >= 0, precision < 10)
  stopifnot(k > r1, precision > 1)

  eff_low <- eff_high <- seq(0, 1, 0.01)
  tmp <- sapply(seq(0, 1, 0.01), function(bla) simon2stage_pval(n1=n1,n2=n-n1,r1=r1,k=k,p0=bla))
  index_low  <- which(tmp >= alpha)[1]
  index_high <- which(tmp >= (1 - alpha))[1] - 1

  if (precision > 2) {
    for (i in 2:precision){
      eff_low  <- seq(eff_low[index_low - 1], eff_low [index_low]     ,10^-i)
      eff_high <- seq(eff_high[index_high]  , eff_high[index_high + 1],10^-i)

      tmp_low  <- sapply(eff_low , function(bla) simon2stage_pval(n1=n1,n2=n-n1,r1=r1,k=k,p0=bla))
      tmp_high <- sapply(eff_high, function(bla) simon2stage_pval(n1=n1,n2=n-n1,r1=r1,k=k,p0=bla))

      index_low  <- which(tmp_low  >= alpha)[1]
      index_high <- which(tmp_high >= (1 - alpha))[1] - 1
    }
  }

  result <- data.frame(CI_low = eff_low[index_low], CI_high = eff_high[index_high])
  result
}
