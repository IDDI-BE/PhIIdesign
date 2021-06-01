
#' @title p-value for a 2-stage PhII single arm Sargent design
#' @description This function calculates the p-value for a 2-stage PhII single arm Sargent design
#' In the case of Sargent's 2-stage design, stage-wise ordering means that outcomes observed in
#' the second stage of the trial are more extreme than outcomes observed in the first stage of the trial.
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
#' sargent2stage_pval(n1=15, n2=10, r1=1, k=6, p0=0.1)

sargent2stage_pval <- function(n1, n2=NULL, r1, k, p0) {
  if (k>r1){
    i <- seq(r1 + 1, n1)
    pval<-sum(dbinom(x=i,size=n1,prob=p0)*(1-pbinom(q=k-i-1,size=n2,prob=p0,lower.tail=T)))
  }
  if (k<=r1 | is.null(n2)){
    pval<-1 - pbinom(q=k-1,size=n1,prob=p0,lower.tail=T)
  }
  return(pval)
}
