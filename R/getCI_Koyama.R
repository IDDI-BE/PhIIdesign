
#' @title Koyama confidence interval.
#' @description
#' Calculates the two sided 1-2*alpha confidence interval based, by inverting the "rejection probability
#' function" (or the "p value function") at every possible value for the null hypothesis between [0,1].
#' A two-sided (1 − 2*alpha) CI interval includes then all values of binomial proportions for
#' which the p value for testing H0, lays within the interval [alpha, 1 − alpha]
#' A stage-wise ordering of the sample space is used:
#' Let m be the stopping stage and x the total number of responses. A trial outcome (m'; x') is at least as
#' extreme (against H0) as the observed trial outcome (m; x), if one of the 3 following conditions is met:
#' (A) m' = m and x'>=x
#' (B) m' = 1, m = 2 and x'>=a
#' (C) m' = 2, m = 1 and x <=r1
#' In the case of Simon’s and Fleming's 2-stage design, stage-wise ordering means that outcomes observed in
#' the second stage of the trial are more extreme than outcomes observed in the first stage of the trial.
#' In the case of Fleming's 2-stage design, stage-wise ordering means that rejection of H0
#' in the first stage is more extreme than rejection in the second stage
#' @param N overall sample size.
#' @param n1 sample size for the first stage.
#' @param r1 critical value for futility decision first stage (Simon's or Sargent's 2-stage)
#' @param a critical value for efficacy decision first stage (Fleming's 2-stage)
#' @param k overall observed responses (must be larger than r1).
#' @param alpha determining the two sided 1-2*alpha confidence interval.
#' @param design 3 choices for 2-stage design: "Simon","Sargent" or "Fleming"
#' @param precision gives the precision (in decimal numbers) to which the confidence interval should be calculated (should be less than 10).
#' @references Koyama T and Chen H (2008): Proper inference from Simon's two-stage designs. Statistics in Medicine, 27(16):3145-3154.
#'     Nhacolo A, Brannath W. Interval and point estimation in adaptive Phase II trials with binary endpoint. Stat Methods Med Res 2019;28:2635-2648
#' @export
#' @examples
#' #Calculate a Simon's two-stage design
#' res <- simon2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.2,eps = 0.005, N_min = 1, N_max = 50)
#' CI <- mapply(a_ = res$N, b_ = res$n1 , c_ = res$r1, d_ = res$r2 + 1, FUN = function(a_, b_, c_, d_)
#' getCI_Koyama(N = a_,  n1 = b_, r1 = c_, k = d_, a=NULL, alpha = 0.05, design="Simon",precision = 3))
#' getCI_Koyama(N = 21,  n1 = 14, r1 = 1, k = 5, a=NULL, alpha = 0.05, design="Simon",precision = 3)

getCI_Koyama<-function (N,n1,r1,a=NULL,k,alpha,design,precision=3){

  stopifnot(class(k)         == "numeric" | class(k) == "integer", k >= 0)
  stopifnot(class(r1)        == "numeric", r1 >= 0)
  stopifnot(class(n1)        == "numeric", n1 >= 0, n1 > r1)
  stopifnot(class(N)         == "numeric", N >= 0, N > n1)
  stopifnot(class(alpha)     == "numeric", alpha > 0, alpha <= 1)
  stopifnot(class(precision) == "numeric", precision >= 0, precision < 10)
  stopifnot(k > r1, precision > 1)

  eff_low <- eff_high <- seq(0, 1, 0.01)

  if (design=="Simon"  ){tmp <- sapply(seq(0, 1, 0.01), function(bla) simon2stage_pval(n1=n1,n2=N-n1,r1=r1,k=k,p0=bla))}
  if (design=="Sargent"){tmp <- sapply(seq(0, 1, 0.01), function(bla) sargent2stage_pval(n1=n1,n2=N-n1,r1=r1,k=k,p0=bla))}
  if (design=="Fleming"){tmp <- sapply(seq(0, 1, 0.01), function(bla) fleming2stage_pval(n1=n1,n2=N-n1,r1=r1,a=a,k=k,p0=bla))}

  index_low  <- which(tmp >= alpha)[1]
  index_high <- which(tmp >= (1 - alpha))[1] - 1

  if (precision > 2) {
    for (i in 2:precision){
      eff_low  <- seq(eff_low[index_low - 1], eff_low [index_low]     ,10^-i)
      eff_high <- seq(eff_high[index_high]  , eff_high[index_high + 1],10^-i)

      if (design=="Simon"){
        tmp_low  <- sapply(eff_low , function(bla) simon2stage_pval  (n1=n1,n2=N-n1,r1=r1,    k=k,p0=bla))
        tmp_high <- sapply(eff_high, function(bla) simon2stage_pval  (n1=n1,n2=N-n1,r1=r1,    k=k,p0=bla))
      }

      if (design=="Sargent"){
        tmp_low  <- sapply(eff_low , function(bla) sargent2stage_pval(n1=n1,n2=N-n1,r1=r1,    k=k,p0=bla))
        tmp_high <- sapply(eff_high, function(bla) sargent2stage_pval(n1=n1,n2=N-n1,r1=r1,    k=k,p0=bla))
      }

      if (design=="Fleming"){
        tmp_low  <- sapply(eff_low , function(bla) fleming2stage_pval(n1=n1,n2=N-n1,r1=r1,a=a,k=k,p0=bla))
        tmp_high <- sapply(eff_high, function(bla) fleming2stage_pval(n1=n1,n2=N-n1,r1=r1,a=a,k=k,p0=bla))
      }

      index_low  <- which(tmp_low  >= alpha)[1]
      index_high <- which(tmp_high >= (1 - alpha))[1] - 1
    }
  }

  result <- data.frame(CI_LL = eff_low[index_low], CI_UL = eff_high[index_high])
  result
}
