
#' @title p-value for a 2-stage PhII single arm Sargent design
#' @description This function calculates the p-value for a 2-stage PhII single arm Fleming's design
#' A stage-wise ordering of the sample space is used:
#' Let m be the stopping stage and x the total number of responses. A trial outcome (m'; x') is at least as
#' extreme (against H0) as the observed trial outcome (m; x), if one of the 3 following conditions is met:
#' (A) m' = m and x'>=x
#' (B) m' = 1, m = 2 and x'>=a
#' (C) m' = 2, m = 1 and x <=r1
#' in other words: rejection of H0 in the first stage is more extreme
#' @param n1 total number of patients in stage1
#' @param n2 total number of patients in stage2
#' @param r1 critical value futility for the first stage
#' @param a critical value efficacy for the first stage
#' @param k overall observed responses (must be larger than r1)
#' @param p0 true success probability under H0
#' @return p-value
#' @references Mander AP, Thompson SG. Two-stage designs optimal under the alternative hypothesis for phase II cancer clinical trials.
#'     Contemporary Clinical Trials 2010;31:572–578
#'     Qin F et al. Optimal, minimax and admissible two-stage design for phase II oncology clinical trials. BMC Medical Research Methodology 2020;20:126
#'     Nhacolo A, Brannath W. Interval and point estimation in adaptive Phase II trials with binary endpoint. Stat Methods Med Res 2019;28:2635-2648
#' @export
#' @examples
#' fleming2stage_pval(n1=21, n2=24, r1=5, a=10,k=17, p0=0.3)

fleming2stage_pval <- function(n1, n2=NULL, r1, a, k, p0){

  if (k>r1){
    nmax <- n1+n2
    b_p0 <- lapply(1:nmax, FUN = function(n) dbinom(0:nmax, size = n, prob = p0))
    B_p0 <- lapply(1:nmax, FUN = function(n) pbinom(0:nmax, size = n, prob = p0, lower.tail = TRUE))

    x1  <- seq(r1 + 1,a-1)

    B_p0_a  <- 1-B_p0[[n1]][1 + a-1]       # P(X1>=a|n1,p0)    so 1-pbinom(a-1   ,n1,p0), indexing with "+1", as B_p0 starts from 0 successes
    b_p0_x1 <-   b_p0[[n1]][1 + x1]        # P(X1=x1|n1,p0)    so   dbinom(x1    ,n1,p0), indexing with "+1", as b_p0 starts from 0 successes
    B_p0_k  <- 1-B_p0[[n2]][1 + (k-1-x1)]  # P(X2>=k-x1|n2,p0) so 1-pbinom(k-1-x1,n2,p0), indexing with "+1", as B_p0 starts from 0 successes
    pval<-B_p0_a+sum(b_p0_x1 * B_p0_k)
  }

  if (k<=r1 | is.null(n2)){
    pval<- 1 - pbinom(q=k-1,size=n1,prob=p0,lower.tail=T)
  }

  return(pval)
}
