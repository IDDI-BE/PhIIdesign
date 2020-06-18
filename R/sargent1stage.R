
#' @title The Sargent 1-stage function
#' @description This function calculates sample sizes of the sargent 1-stage design.
#' @details if x<=r --> stop futility
#' @details if x>=s --> efficacy
#' @param p0 uninteresting response (null hypothesis H0)
#' @param pa interesting response (alternative hypothesis Ha)
#' @param alpha P(reject H0|H0)
#' @param beta P(reject Ha|Ha)
#' @param eta P(reject Ha|H0)
#' @param pi P(reject H0|Ha)
#' @param eps tolerance (actual alpha<=alpha+eps; actual beta<=beta+eps; actual eta>=eta-eps; actual p>=pi-eps); default value = 0.005
#' @param N_min minimum sample size value for grid search
#' @param N_max maximum sample size value for grid search
#' @examples
#' sargent1stage(p0=0.5,pa=0.65,alpha=0.1,beta=0.1,eta=0.8,pi=0.8,eps = 0.005,N_min=0,N_max=100)
#' @export
#' @import data.table
#' @importFrom stats aggregate dbinom

# p0=0.2;pa=0.35;alpha=0.1;beta=0.1;eta=0.8;pi=0.8;eps = 0.005;N_min=35;N_max=50

sargent1stage<-function(p0,pa,alpha,beta,eta,pi,eps = 0.005,N_min,N_max){

  if (pa<p0) {stop('p0 should be smaller than pa')}

  # Define variables as a NULL value (to avoid 'notes' in devtools package check)

  N<-NULL

  # Get for all N's, all possible r's (0:N-2) and select where beta_temp=P(X<=r|Ha)<=beta+eps
  #------------------------------------------------------------------------------------------

  res0           <- data.frame(do.call("rbind", mapply(function(a) cbind(N=a,r=0:(a-2)),a=(N_min:N_max),SIMPLIFY=F)))
  res0$beta_temp <- pbinom(q=res0$r,size=res0$N,prob=pa,lower.tail=T)
  res1           <- res0[res0$beta_temp<=(beta +eps),]
  names(res1)    <- c("N","r","beta_temp")

  # Calculate eta for all scenarios (only r needed)
  #------------------------------------------------

  res1$eta_temp  <- pbinom(q=res1$r,size=res1$N,prob=p0,lower.tail=T)
  res1$diff_eta  <- eta - res1$eta_temp
  res2           <- res1[res1$diff_eta<=eps,]

  # Get for selected N's, r's: s's ((r+2):n), and calculate pi
  #----------------------------------------------------------

  res3          <- data.frame(do.call("rbind", mapply(function(a,b,c,d,e,f,g) cbind(N=a,r=b,beta_temp=c,eta_temp=d,s=(b+2):a),
                                                                                    a=res2$N,b=res2$r,c=res2$beta_temp,d=res2$eta_temp)))
  # Calculate pi all scenarios (s needed)
  #--------------------------------------

  res3$pi_temp  <- 1-pbinom(q=res3$s-1,size=res3$N,prob=pa,lower.tail=T)
  res3$diff_pi  <- pi - res3$pi_temp
  res4          <- res3[res3$diff_pi<=eps,]

  # Calculate alpha for all scenarios (s needed)
  #--------------------------------------------

  res4$alpha_temp<- 1-pbinom(q=res4$s-1,size=res4$N,prob=p0,lower.tail=T)
  res4$diff_alpha<- res4$alpha_temp - alpha
  res5           <- res4[res4$diff_alpha<=eps,]

  res5$alpha <- alpha
  res5$beta  <- beta
  res5$eta   <- eta
  res5$pi    <- pi

  res5 <-data.table::as.data.table(res5)
  res5 <- res5[,N_min := min(N) ]
  res  <-  res5[which(N==N_min),]

  names(res)[names(res)=="alpha"]<-"alpha_param"
  names(res)[names(res)=="beta" ]<-"beta_param"
  names(res)[names(res)=="eta"  ]<-"eta_param"
  names(res)[names(res)=="pi"   ]<-"pi_param"

  names(res)[names(res)=="alpha_temp"]<-"alpha"
  names(res)[names(res)=="beta_temp" ]<-"beta"
  names(res)[names(res)=="eta_temp"  ]<-"eta"
  names(res)[names(res)=="pi_temp"   ]<-"pi"

  res<- res[,-c("diff_pi","diff_alpha","N_min")]
  res<- cbind(res[,c("r","s","N","alpha","beta","eta","pi")],p0=p0,pa=pa,res[,c("alpha_param","beta_param","eta_param","pi_param")])

  return(res)

}

# TEST (Sargent DJ, Goldberg RM. A Three-Outcome Design for Phase II Clinical Trials. Controlled Clinical Trials 22:117-125
#--------------------------------------------------------------------------------------------------------------------------
#test_0<- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,2),pa=rep((a+0.15),2),alpha=c(0.1,0.05),beta=c(0.1,0.1),eta=c(0.8,0.8),pi=c(0.8,0.8)),a=c(0.05,0.1,0.2,0.3,0.4,0.5),SIMPLIFY=F)))
#
#for (i in 1:dim(test_0)[1]){
#  res<- PhIIdesign::sargent1stage(p0=test_0[i,]$p0,pa=test_0[i,]$pa,alpha=test_0[i,]$alpha,beta=test_0[i,]$beta,eta=test_0[i,]$eta,pi=test_0[i,]$pi,eps=0.005,N_min=20,N_max=70)
#
#  if (i==1) {test_list      <-list(res)}   # Create list with first dataset
#  if (i!=1) {test_list[[i]] <-res }        # Next iterations: append dataset
#
#}
#test<- data.frame(do.call("rbind",test_list))
#test<-test[,c("p0","pa","alpha_param","beta_param","eta_param","pi_param","alpha","beta","eta","pi","r","s","N")]
