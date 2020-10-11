
#' @title P(reject Ha) in a 2-stage PhII single arm Simon design
#' @description This function calculates the P(reject alternative hypothesis|true success probability)
#' @param n1 total number of patients in stage1
#' @param n2 total number of patients in stage2
#' @param r1 critical value for the first stage
#' @param r critical value for the second stage
#' @param p true success probability
#' @return Probability
#' @details
#' if x1<=r1 --> stop futility \cr
#' if (x1+x2)<=r --> futility \cr
#' with x1 the number of successes in the first stage and x2 the number of successes in the second stage
#' @references Simon R. Optimal two-Stage Designs for Phase II Clinical Trials. \cr
#' Control Clin Trials. 1989;10:1-10
#' @export
#' @examples
#' simon2stage_prob_reject_Ha(n1=10,n2=20,r1=0,r=4,p=0.25)

simon2stage_prob_reject_Ha <- function(n1, n2, r1, r, p) {
  i <- seq(r1 + 1, min(n1, r))
  pbinom(q = r1, size = n1, prob = p, lower.tail = TRUE) +
    sum(dbinom(i, size = n1, prob = p) * pbinom(q = r - i, size = n2, prob = p, lower.tail = TRUE))
}
