

probsimon<-function(n1,n2,r1,r,p){
  i<-seq(r1+1,min(n1,r))
  pbinom(q=r1,size=n1,prob=p,lower.tail=T)+sum(dbinom(i,size=n1,prob=p)*pbinom(q=r-i,size=n2,prob=p,lower.tail=T))
}


#' @title The Simon 2-stage function
#' @description This function calculates sample sizes of the sargent 2-stage design.
#' @details n1= n stage1
#' @details n2= n stage2
#' @details total number of patients=n1+n2
#' @details if x1<=r1 --> stop futility
#' @details if (x1+x2)<=r --> futility
#' @details if (x1+x2)> s --> efficacy
#' @details EN.p0 = expected sample size under H0
#' @details PET.p0 = probability of terminating the trial at the end of the first stage under H0
#' @details Confidence interval according to Koyama T, Chen H. Proper inference from simon’s two-stage designs. Stat Med. 2008; 27:3145–154;
#'
#' @param p0 uninteresting response (null hypothesis H0)
#' @param pa interesting response (alternative hypothesis Ha)
#' @param alpha P(reject H0|H0)
#' @param beta P(reject Ha|Ha)
#' @param eps tolerance (actual alpha<=alpha+eps; actual beta<=beta+eps; actual eta>=eta-eps; actual p>=pi-eps); default value = 0.005
#' @param N_min minimum sample size value for grid search
#' @param N_max maximum sample size value for grid search
#' @examples
#' simon2stage(p0=0.1,pa=0.3,alpha=0.05,beta=0.2,eps = 0.005,N_min=0,N_max=50)
#' @export
#' @import data.table
#' @import OneArmPhaseTwoStudy
#' @importFrom grDevices chull
#' @importFrom graphics plot points
#' @importFrom stats aggregate dbinom

simon2stage<-function(p0,pa,alpha,beta,eps = 0.005,N_min,N_max){

  if (pa<p0) {stop('p0 should be smaller than pa')}

  # Define variables as a NULL value (to avoid 'notes' in devtools package check)

  EN.p0<-EN.p0_N_min<-EN.p0_min<-MIN<-N<-OPT<-NULL

  #----------------------------------------------------------------------------------------------------------#
  # Get all possible scenarios for N, n1, r1, r2 and n2                                                      #
  # Note that this is the possible range of values: 0 <=r1<n1; r1+1<=r; r<=n1+n2                             #
  #----------------------------------------------------------------------------------------------------------#

  # Get all N"s for which there is a max r (1:N) for which P(X<=r|Ha)<=beta+eps, and select that max r
  # Note that this is not taking into account first stage. However, the probability of rejecting Ha
  # should be lower in the second stage, compared to a 1-stage design, as some cases already rejected at first stage,
  # meaning that an rmax, calculated using a 1-stage design, is a good maximum
  #------------------------------------------------------------------------------------------------------------

  res0_t           <- data.frame(do.call("rbind", mapply(function(a) cbind(N=a,rtemp=1:a),a=(N_min:N_max),SIMPLIFY=F)))
  res0_t$betamax   <- pbinom(q=res0_t$rtemp,size=res0_t$N,prob=pa,lower.tail=T)
  res0             <- aggregate(res0_t[!is.na(res0_t$betamax) & res0_t$betamax<=(beta +eps),]$rtemp, by = list(res0_t[!is.na(res0_t$betamax) & res0_t$betamax<=(beta +eps),]$N), max)    # Select all possible N"s: there needs to be
                                                                                                                                                                                         # at least one r with P(X<=r|Ha)<=beta+eps
  names(res0)      <- c("N","rmax")

  # Get for selected N's all possible n1's (1:N-1) + create r1max with 0<=r1<=r-1 and 0<=r1<n1 (note: r1<n1, otherwise first stage makes
  # no sense: futility even if all outcomes a success)
  #-------------------------------------------------------------------------------------------

  res1        <- data.frame(do.call("rbind", mapply(function(a,b) cbind(N=a,n1=(1:(a-1)),rmax=b),a=res0$N,b=res0$rmax,SIMPLIFY=F)))

  res1$r1max  <- pmin(res1$rmax-1,res1$n1-1)

  # Get for selected N's and n1's, r1's (0:r1max, where P(X_1<=r_1|Ha)<=beta+eps)
  #------------------------------------------------------------------------------

  res2_t      <- data.frame(do.call("rbind", mapply(function(a,b,c,d) cbind(N=a,n1=b,rmax=c,r1max=d,r1=(0:d)),
                                                                            a=res1$N,b=res1$n1,c=res1$rmax,d=res1$r1max,SIMPLIFY=F)))

  res2_t$beta <- pbinom(q=res2_t$r1,size=res2_t$n1,prob=pa,lower.tail=T)
  res2        <- res2_t[res2_t$beta<=(beta +eps),]

  # Get for selected N"s, n1"s and r1"s: r2"s ((r1+1):rmax)
  #----------------------------------------------------------

  res3        <- data.frame(do.call("rbind", mapply(function(a,b,c,d) cbind(N=a,n1=b,r1=d,r2=((d+1):c)),
                                                                            a=res2$N,b=res2$n1,c=res2$rmax,d=res2$r1)))
  res3$n2     <- res3$N-res3$n1


  # Calculate beta and alpha
  #-------------------------

  res3$beta_temp   <- mapply(function (a,b,c,d) probsimon(n1=a,n2=b,r1=c,r=d,p=pa),a=res3$n1,b=res3$n2,c=res3$r1,d=res3$r2)
  res3$diff_beta   <- res3$beta_temp - beta
  res4             <- res3[res3$diff_beta<=eps,]

  res4$alpha_temp  <- 1-mapply(function (a,b,c,d) probsimon(n1=a,n2=b,r1=c,r=d,p=p0),a=res4$n1,b=res4$n2,c=res4$r1,d=res4$r2)
  res4$diff_alpha  <- res4$alpha_temp - alpha
  res5             <- res4[res4$diff_alpha<=eps,]

  res5$alpha <- alpha
  res5$beta  <- beta

  res5$PET.p0 <- pbinom(q=res5$r1,size=res5$n1,prob=p0,lower.tail=T)
  res5$EN.p0  <- res5$N-((res5$N-res5$n1)*res5$PET.p0)

  res5<-data.table::as.data.table(res5)

  res5 <- res5[,N_min := min(N) ]
  res5 <- res5[,EN.p0_min := min(EN.p0)]
  res5 <- res5[,EN.p0_N_min := min(EN.p0),by=N]
  res6 <- res5[EN.p0==EN.p0_N_min]

  res6$OPT<-res6$MIN<-res6$ADMISS<-c("")
  res6[which(EN.p0==EN.p0_min),]$OPT  <-"Optimal"
  res6[which(N==N_min),]$MIN          <-"Minimax" # Note: if multiple designs that meet the criteria for minimax:choose
  # design with minimal expected sample size under H0: "Optimal minimax design"

  # Get admissible designs
  y <- data.frame(res6[,c("N","EN.p0")])
  con.ind <- chull(y)[chull((y)) == cummin(chull((y)))]

  #chull_result<-data.frame(print(multichull::CHull(y,bound = "lower")))
  #con.ind <- as.numeric(rownames(y[y$N %in% chull_result$complexity,]))

  #plot(y$N, y$EN.p0)
  #lines(y[con.ind,]$N, y[con.ind,]$EN.p0)

  res<-res6[N>=min(res6[MIN=="Minimax","N"]) & N<=max(res6[OPT=="Optimal","N"]),]
  res[which((rownames(res) %in% c(con.ind)) & (N > N_min) &  (EN.p0 > EN.p0_min) & (N<res[OPT=="Optimal",]$N)),]$ADMISS<- "Admissible"

  names(res)[names(res)=="alpha"]<-"alpha_param"
  names(res)[names(res)=="beta" ]<-"beta_param"
  names(res)[names(res)=="eta"  ]<-"eta_param"
  names(res)[names(res)=="pi"   ]<-"pi_param"

  names(res)[names(res)=="alpha_temp"]<-"alpha"
  names(res)[names(res)=="beta_temp" ]<-"beta"
  names(res)[names(res)=="eta_temp"  ]<-"eta"
  names(res)[names(res)=="pi_temp"   ]<-"pi"

  # Plot
  #-----

  res.df<-data.frame(res)
  plot(res.df[,"N"],res.df[,"EN.p0"], type = "l",
       xlab = "Maximum Sample Size N",
       ylab = expression(paste("E( N | ",p[0], " )")),
       main = "Two-stage Designs")
  points(res.df[res.df$MIN   =="Minimax"   , "N"],res.df[res.df$MIN   =="Minimax"   , "EN.p0"],pch = "M")
  points(res.df[res.df$OPT   =="Optimal"   , "N"],res.df[res.df$OPT   =="Optimal"   , "EN.p0"],pch = "O")
  points(res.df[res.df$ADMISS=="Admissible", "N"],res.df[res.df$ADMISS=="Admissible", "EN.p0"],pch = "A")

  # Calculate 1-2*alpha confidence interval, based on Koyama, Statistics in Medicine 2008

  res$eff<-paste0(res$r2+1,"/",res$N," (",100*round((res$r2+1)/res$N,3),"%)")

  CI<-mapply(function(a,b,c,d) OneArmPhaseTwoStudy::get_CI(k=a,r1=b,n1=c,n=d,alpha=alpha,precision=3),
                                                           a=res$r2+1,b=res$r1,c=res$n1,d=res$N)
  res$CI_low  <-100*unlist(CI[rownames(CI)=="CI_low" ,])
  res$CI_high <-100*unlist(CI[rownames(CI)=="CI_high",])

  res<- res[,-c("diff_alpha","N_min","EN.p0_min","EN.p0_N_min")]
  res<- cbind(res[,c("r1","n1","r2","n2","N","eff","CI_low","CI_high","EN.p0","PET.p0","MIN","OPT","ADMISS","alpha","beta")],p0=p0,pa=pa,res[,c("alpha_param","beta_param")])

  names(res)[names(res)=="CI_low" ]<- paste0(100-2*100*alpha,"%CI_low")
  names(res)[names(res)=="CI_high"]<- paste0(100-2*100*alpha,"%CI_high")

  res<<-res
  return(res)

}

# TEST (Simon R. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 1989;10:1-10 : Table 1 and 2
#-------------------------------------------------------------------------------------------------------------------------------
#test_1_0  <- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,3),pa=rep((a+0.2),3),alpha=c(0.1,0.05,0.05),beta=c(0.1,0.2,0.1)),a=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7),SIMPLIFY=F)))
#for(i in 1:dim(test_1_0)[1]){
#  nmax   <- PhIIdesign::fleming1stage(p0=test_1_0[i,]$p0,pa=test_1_0[i,]$pa,alpha=test_1_0[i,]$alpha,beta=test_1_0[i,]$beta,eps=0)$n+15
#  res    <- simon2stage  (p0=test_1_0[i,]$p0,pa=test_1_0[i,]$pa,alpha=test_1_0[i,]$alpha,beta=test_1_0[i,]$beta,eps=0,N_min=10,N_max=nmax)
#  res_dt <- cbind(res[OPT=="Optimal",c("r1","n1","r2","N","EN.p0","PET.p0")],res[MIN=="Minimax",c("r1","n1","r2","N","EN.p0","PET.p0")])
#  if (i==1) {test1_list      <- list(res_dt)}
#  if (i!=1) {test1_list[[i]] <- res_dt }
#}
#test1             <- cbind(test_1_0,data.frame(do.call("rbind",test1_list)))
#
#test_2_0  <- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,3),pa=rep((a+0.15),3),alpha=c(0.1,0.05,0.05),beta=c(0.1,0.2,0.1)),a=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),SIMPLIFY=F)))
#for(i in 1:dim(test_2_0)[1]){
#  nmax   <- PhIIdesign::fleming1stage(p0=test_2_0[i,]$p0,pa=test_2_0[i,]$pa,alpha=test_2_0[i,]$alpha,beta=test_2_0[i,]$beta,eps=0)$n+20
#  res    <- PhIIdesign::simon2stage  (p0=test_2_0[i,]$p0,pa=test_2_0[i,]$pa,alpha=test_2_0[i,]$alpha,beta=test_2_0[i,]$beta,eps=0,N_min=25,N_max=nmax)
#  res_dt <- cbind(res[OPT=="Optimal",c("r1","n1","r2","N","EN.p0","PET.p0")],res[MIN=="Minimax",c("r1","n1","r2","N","EN.p0","PET.p0")])
#  if (i==1) {test2_list      <- list(res_dt)}
#  if (i!=1) {test2_list[[i]] <- res_dt }
#}
#test2             <- cbind(test_2_0,data.frame(do.call("rbind",test2_list)))



