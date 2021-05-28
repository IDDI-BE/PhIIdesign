
#' @title p-value for a 2-stage PhII single arm Simon design
#' @description This function calculates the p-value for a 2-stage PhII single arm Simon design
#' @param n1 total number of patients in stage1
#' @param n2 total number of patients in stage2
#' @param r1 critical value for the first stage
#' @param k overall observed responses (must be larger than r1)
#' @param p0 true success probability under H0
#' @return p-value
#' @references Simon R. Optimal two-Stage Designs for Phase II Clinical Trials. \cr
#' Control Clin Trials. 1989;10:1-10; Koyama T and Chen H (2008): Proper inference from Simon's two-stage designs. \cr
#' Statistics in Medicine, 27(16):3145-3154.
#' @export
#' @examples
#' simon2stage_pval(n1=15, n2=10, r1=1, k=6, p0=0.1)
#' #Note: same result with:
#' 1-simon2stage_prob_reject_Ha(n1=15, n2=10, r1=1, r=5, p=0.1)

simon2stage_pval <- function(n1, n2, r1, k, p0) {
  if (k>r1){
    i <- seq(r1 + 1, n1)
    sum(dbinom(x=i,size=n1,prob=p0)*(1-pbinom(q=k-i-1,size=n2,prob=p0,lower.tail=T)))
  }
  else{
    1 - pbinom(q=k-1,size=n1,prob=p0,lower.tail=T)
  }
}
