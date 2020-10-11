

#' @title P(reject H0) in a 2-stage PhII single arm Sargent design
#' @description This function calculates the P(reject null hypothesis|true success probability)
#' @param n1 total number of patients in stage1
#' @param n2 total number of patients in stage2
#' @param r1 critical value for the first stage
#' @param s critical value for the second stage
#' @param p true success probability
#' @return Probability
#' @details
#' if x1<=r1 --> stop futility \cr
#' if (x1+x2)>=s --> success \cr
#' with x1 the number of successes in the first stage and x2 the number of successes in the second stage
#' @references Sargent DJ, Chan V, Goldberg RM. A three-outcome design for phase II clinical \cr
#' trials. Control Clin Trials. 2001;22(2):117-125. doi:10.1016/s0197-2456(00)00115-x
#' @export
#' @examples
#' sargent2stage_prob_reject_H0(n1=10,n2=20,r1=0,s=4,p=0.25)

sargent2stage_prob_reject_H0<-function(n1,n2,r1,s,p){
  i<-seq(r1+1,min(n1,s))
  sum(pbinom(q=s-i-1,size=n2,prob=p,lower.tail=F)*dbinom(i,size=n1,prob=p))+as.numeric(n1>s)*pbinom(q=s,size=n1,prob=p,lower.tail=F)
}
