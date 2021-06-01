

#' @title P(reject Ha) in a 2-stage PhII single arm Fleming design
#' @description This function calculates the P(reject alternative hypothesis)
#' @param n1 total number of patients in stage1
#' @param n2 total number of patients in stage2
#' @param r1 critical value futility for the first stage
#' @param a critical value efficacy for the first stage
#' @param r2 critical value futility/efficacy for the second stage
#' @param p0 true success probability under H0
#' @param pa true success probability under Ha
#' @details
#' if x1<=r1 --> stop futility \cr
#' if x1>=a --> stop efficacy \cr
#' if (x1+x2)<=r2 --> futility \cr
#' with x1 the number of successes in the first stage and x2 the number of successes in the second stage
#' @references Mander AP, Thompson SG. Two-stage designs optimal under the alternative hypothesis for phase II cancer clinical trials.
#'     Contemporary Clinical Trials 2010;31:572–578
#'     Qin F et al. Optimal, minimax and admissible two-stage design for phase II oncology clinical trials. BMC Medical Research Methodology 2020;20:126
#' @export
#' @examples
#' fleming2st_reject_Ha(n1=13, n2=7, r1=0, a=3, r2=2, p0=0.05, pa=0.25)


fleming2st_reject_Ha <- function(n1, n2, r1, a, r2, p0, pa){
  ## r2 is the total number of seen cases on both n1 and n2 together
  ## Calculate alpha/beta for all values up to r2

  nmax <- n1+n2
  b_p0 <- lapply(1:nmax, FUN = function(n) dbinom(0:nmax, size = n, prob = p0))
  b_pa <- lapply(1:nmax, FUN = function(n) dbinom(0:nmax, size = n, prob = pa))
  B_p0 <- lapply(1:nmax, FUN = function(n) pbinom(0:nmax, size = n, prob = p0, lower.tail = TRUE))
  B_pa <- lapply(1:nmax, FUN = function(n) pbinom(0:nmax, size = n, prob = pa, lower.tail = TRUE))

  x1  <- seq(r1 + 1,a-1)

  # calculate eta_temp=P(not rejecting H0) under H0 for each (r1,x1)
  #-----------------------------------------------------------------

  B_p0_r1    <- B_p0[[n1]][1 + r1]       # P(X1<=r1|n1,p0)   so pbinom(r1 ,n1,p0), indexing with "+1", as B_p0 starts from 0 successes
  b_p0_x1    <- b_p0[[n1]][1 + x1]       # P(X1=x1|n1,p0)    so dbinom(x1  ,n1,p0), indexing with "+1", as b_p0 starts from 0 successes
  B_p0_r2    <- B_p0[[n2]][1 + (r2-x1)]  # P(X2<=r2-x1|n2,p0) so pbinom(r2-x1,n2,p0), indexing with "+1", as B_p0 starts from 0 successes
  eta_temp   <- B_p0_r1 + sum(b_p0_x1 * B_p0_r2)

  # calculate beta_temp=P(not rejecting H0) under Ha for each (r1,x1)
  #------------------------------------------------------------------

  B_pa_r1    <- B_pa[[n1]][1 + r1]
  b_pa_x1    <- b_pa[[n1]][1 + x1]
  B_pa_r2    <- B_pa[[n2]][1 + (r2-x1)]
  beta_temp  <- B_pa_r1 + sum(b_pa_x1 * B_pa_r2)

  PET.p0 <- pbinom(q = r1, size = n1, prob = p0, lower.tail = TRUE) + pbinom(q = a-1, size = n1, prob = p0, lower.tail = FALSE)
  PET.pa <- pbinom(q = r1, size = n1, prob = pa, lower.tail = TRUE) + pbinom(q = a-1, size = n1, prob = pa, lower.tail = FALSE)

  EN.p0  <- (n1+n2) - (n2*PET.p0)
  EN.pa  <- (n1+n2) - (n2*PET.pa)

  return(list(N=n1+n2,n1=n1,r1=r1,a=a,n2=n2,r2=r2,alpha_temp=1-eta_temp,beta_temp=beta_temp,EN.p0=EN.p0,EN.pa=EN.pa,PET.p0=PET.p0,PET.pa=PET.pa))
}
