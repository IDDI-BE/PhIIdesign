
#' @title P(reject Ha) in a 2-stage PhII single arm Sargent design
#' @description This function calculates the P(reject alternative hypothesis|true success probability)
#' @param n1 total number of patients in stage1
#' @param n2 total number of patients in stage2
#' @param r1 critical value for the first stage
#' @param r2 critical value for the second stage
#' @param p true success probability
#' @return Probability
#' @details
#' if x1<=r1 --> stop futility \cr
#' if (x1+x2)<=r --> futility \cr
#' with x1 the number of successes in the first stage and x2 the number of successes in the second stage
#' @references Sargent DJ, Chan V, Goldberg RM. A three-outcome design for phase II clinical trials. Control Clin Trials. 2001;22(2):117-125. doi:10.1016/s0197-2456(00)00115-x
#' @export
#' @examples
#' sargent2stage_prob_reject_Ha(n1=10,n2=20,r1=0,r2=4,p=0.25)

sargent2stage_prob_reject_Ha<-function(n1,n2,r1,r2,p){
  i<-seq(r1+1,min(n1,r2))
  pbinom(q=r1,size=n1,prob=p,lower.tail=T)+sum(pbinom(q=r2-i,size=n2,prob=p,lower.tail=T)*dbinom(i,size=n1,prob=p))
}

