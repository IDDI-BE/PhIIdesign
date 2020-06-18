
#-----------------------------------------------------------------------------------------#
# probability distribution function for difference between two binomial variables Z=X1-X2 #
# (allowing for unequal sample sizes in both groups)                                      #
#-----------------------------------------------------------------------------------------#

bin_dif_pdf<-function(z,n1,n2,p1,p2,type){
  if (type=="exact"){
    i=(0:max(n1,n2))
    prob<- sum(dbinom(x=i+z,size=n1,prob=p1)*dbinom(x=i,size=n2,prob=p2))
  }

  if (type=="normal"){
    prob<-pnorm(q=z+0.5,mean=(n1*p1-n2*p2),sd=(n1*p1*(1-p1)+n2*p2*(1-p2))**0.5,lower.tail=T) -
          pnorm(q=z-0.5,mean=(n1*p1-n2*p2),sd=(n1*p1*(1-p1)+n2*p2*(1-p2))**0.5,lower.tail=T)
  }
  return(prob)
}


#----------------------------------------------------------------------------------------#
# Cumulative distribution function for difference between two binomial variables Z=X1-X2 #
# (allowing for unequal sample sizes in both groups)                                     #
#----------------------------------------------------------------------------------------#

bin_dif_cdf<-function(z,n1,n2,p1,p2,type){  # z=cutoff difference between two binomial variables; Z=X1-X2

  if (type=="exact"){
    prob=0

    for (diff in -n2:z){                    # Z=X1-X2, so -n2 is minimum
      i=(0:max(n1,n2))
      prob<-prob+sum(dbinom(x=i+diff,size=n1,prob=p1)*dbinom(x=i,size=n2,prob=p2))
    }
  }

  if (type=="normal"){
    prob<-pnorm(q=z+0.5,mean=(n1*p1-n2*p2),sd=(n1*p1*(1-p1)+n2*p2*(1-p2))**0.5,lower.tail=T)
  }

  return(prob)
}


#-----------------------------------------------------#
# two error probabilities functions (Equations p3530) #
#-----------------------------------------------------#

probho_R<-function(n1_E_,n1_C_,n2_E_,n2_C_,r1_,r_,p_E_,p_C_,...){  # Reject Ha / Accept H0

  # part I of formula

  part1    <- mapply(function(a,b,c) bin_dif_cdf(z=c,n1=a,n2=b,p1=p_E_,p2=p_C_,...),
                     a=n1_E_,b=n1_C_,c=r1_,
                     SIMPLIFY=T)

  # part II of formula
  df      <- data.frame(do.call("rbind", mapply(function(a,b,c,d,e,f,g) cbind(n1_E=a,n1_C=b,n2_E=c,n2_C=d,r1=e,r=f,i=(e+1):g),
                                                a=n1_E_,b=n1_C_,c=n2_E_,d=n2_C_,e=r1_,f=r_,g=pmin(n1_E_,n2_C_+r_),
                                                SIMPLIFY=F)))

  df$res  <- mapply(function(a,b,c,d,e,f) bin_dif_cdf(z=e-f,n1=c,n2=d,p1=p_E_,p2=p_C_,...)*bin_dif_pdf(z=f,n1=a,n2=b,p1=p_E_,p2=p_C_,...),
                    a=df$n1_E,b=df$n1_C,c=df$n2_E,d=df$n2_C,e=df$r,f=df$i,
                    SIMPLIFY=T)

  part2   <- aggregate(df$res, by = list(df$n1_E,df$n1_C,df$n2_E,df$n2_C,df$r1,df$r), FUN=sum)

  probho<-part1+part2$x
  return(probho)
}


probha_R<-function(n1_E_,n1_C_,n2_E_,n2_C_,r1_,s_,p_E_,p_C_,...){  # Reject H0 / Accept Ha

  # part I of formula
  df      <- data.frame(do.call("rbind", mapply(function(a,b,c,d,e,f,g) cbind(n1_E=a,n1_C=b,n2_E=c,n2_C=d,r1=e,s=f,i=(e+1):g),
                                                a=n1_E_,b=n1_C_,c=n2_E_,d=n2_C_,e=r1_,f=s_,g=pmin(n1_E_,n2_C_+s_),
                                                SIMPLIFY=F)))

  df$res  <- mapply(function(a,b,c,d,e,f) (1-bin_dif_cdf(z=e-f-1,n1=c,n2=d,p1=p_E_,p2=p_C_,...))*bin_dif_pdf(z=f,n1=a,n2=b,p1=p_E_,p2=p_C_,...),
                    a=df$n1_E,b=df$n1_C,c=df$n2_E,d=df$n2_C,e=df$s,f=df$i,
                    SIMPLIFY=T)

  part1   <- aggregate(df$res, by = list(df$n1_E,df$n1_C,df$n2_E,df$n2_C,df$r1,df$s), FUN=sum)

  # part II of formula

  part2    <- mapply(function(a,b,c,d)  d*(1-bin_dif_cdf(z=c,n1=a,n2=b,p1=p_E_,p2=p_C_,...)),  # 'c' rather than 'c-1'
                     a=n1_E_,b=n1_C_,c=n2_C_+s_,d=n1_C_>(n2_C_+s_),
                     SIMPLIFY=T)

  probha<-part1$x+part2
  return(probha)
}

#' @title The sargent 2-stage function for two randomized arms

#' @description This function calculates sample sizes of the sargent 2-stage design (randomized)
#' @description including first stage and second stage cutoff values for futility (r1 and r2) and efficacy (s)

#' @param p0 uninteresting response (null hypothesis H0)
#' @param pa interesting response (alternative hypothesis Ha)
#' @param alpha P(reject H0|H0)
#' @param beta P(reject Ha|Ha)
#' @param eta P(reject Ha|H0)
#' @param pi P(reject H0|Ha)
#' @param eps tolerance (actual alpha<=alpha+eps; actual beta<=beta+eps; actual eta>=eta-eps; actual pi>=pi-eps); default value = 0.005
#' @param N_min minimum sample size value for grid search
#' @param N_max maximum sample size value for grid search
#' @param alloc allocation ratio (e.g. 2 for 2:1)
#' @param ... refers to type="exact" or "normal" for pdf and cdf difference of two binomial variables

#' @details E_1= number of successes experimental arm, stage 1
#' @details C_1= number of successes control arm, stage 1
#' @details E_2= number of successes experimental arm, stage 2
#' @details C_2= number of successes control arm, stage 2
#' @details with E_i~Bin(n_i,P) and C_i~Bin(n_i,P), i=1,2
#' @details if (E_1-C_1)<=r1 -> stop futility
#' @details if (E_2+E_1) - (C_2+C_1)<=r2 -> futility
#' @details if (E_2+E_1) - (C_2+C_1)>=s -> efficacy
#' @details Note that all sample sizes can be entered, but only calculations for sample sizes, without any decimals (taking into account allocation ratio)
#' @details Variable of interest: Z=E-C, so always in this direction
#' @details alloc >=1 (more patients in the experimental arm)

#' @return dataframe with selected scenarios, with following columns:
#' \describe{
#'     \item{r1}{if (E_1-C_1)<=r1 -> stop for futility}
#'     \item{n1_E}{Number of patients for experimental arm in stage 1}
#'     \item{n1_C}{Number of patients for control arm in stage 1}
#'     \item{n1}{Total number of patients in stage 1}
#'     \item{r2}{if (E_2+E_1)-(C_2+C_1)<=r1 -> futility}
#'     \item{s}{if (E_2+E_1)-(C_2+C_1)>=1 -> efficacy}
#'     \item{n2_E}{Number of patients for experimental arm in stage 1}
#'     \item{n2_C}{Number of patients for control arm in stage 1}
#'     \item{n2}{Total number of patients in stage 2}
#'     \item{N}{Total number of patients}
#'     \item{N_E}{Total number of patients in experimental arm}
#'     \item{N_C}{Total number of patients in control arm}
#'     \item{EN.p0}{expected sample size, under the null hypothesis}
#'     \item{PET.p0}{robability of terminating the trial at the end of the first stage, under the null hypothesis}
#'     \item{MIN}{Selected 'minimax' design}
#'     \item{OPT}{Selected 'optimal' design}
#'     \item{ADMISS}{Selected 'admissible' design}
#'     \item{alpha}{calculated alpha for selected design}
#'     \item{beta}{calculated beta for selected design}
#'     \item{eta}{calculated eta for selected design}
#'     \item{pi}{calculated pi for selected design}
#'     \item{lambda}{calculated lambda for selected design}
#'     \item{delta}{calculated delta for selected design}
#'     \item{alloc}{allocation rate}
#'     \item{p0}{p0 parameter for design}
#'     \item{pa}{pa parameter for design}
#'     \item{alpha_param}{alpha parameter for design}
#'     \item{beta_param}{beta parameter for design}
#'     \item{eta_param}{eta parameter for design}
#'     \item{pi_param}{pi parameter for design}
#' }

#' @examples
#' \dontrun{
#' sargent2stage_R(p0=0.1,pa=0.3,alpha=0.05,beta=0.15,eta=0.8,pi=0.8,eps = 0.005,
#'     N_min=88,N_max=90,alloc=1,type="normal")
#' }

#' @export
#' @import data.table
#' @importFrom grDevices chull
#' @importFrom graphics plot points
#' @importFrom stats aggregate dbinom pnorm


#------------------#
# Actual function  #
#------------------#
#p0=0.1;pa=0.3;alpha=0.05;beta=0.1;eta=0.85;pi=0.8;eps = 0.005;N_min=100;N_max=104;alloc=1
sargent2stage_R <-function(p0,pa,alpha,beta,eta,pi,eps = 0.005,N_min,N_max,alloc,...) {     # eps=tolerance limit

  if (pa<p0) {stop('p0 should be smaller than pa')}

  # Define variables as a NULL value (to avoid 'notes' in devtools package check)

  EN.p0<-EN.p0_N_min<-EN.p0_min<-MIN<-N<-OPT<-NULL

  p_C           <- p0
  p_E           <- pa
  p_M           <-(p0+pa)/2
  SS_df         <- data.frame(N_E=(N_min:N_max)*(alloc/(alloc+1)),N_C=(N_min:N_max)*(1/(alloc+1)))
  SS_df_rounded <- SS_df[SS_df$N_E==round(SS_df$N_E) & SS_df$N_C==round(SS_df$N_C),]
  N_E           <- SS_df_rounded$N_E
  N_C           <- SS_df_rounded$N_C
  SS            <- N_E+N_C

  #----------------------------------------------------------------------------------------------------------------------------#
  # Get all possible scenarios for N=n1_E+n1_C, r1, r2, n2_E and n2_C, with N=N_E+N_C, with N_E=n1_E+n2_E, and n_C=n1_C+n2_C #
  # Note that this is the possible range of values: -n1_C<=r1<n1_E; -N_C<=rmax<=N_E-2; r+2<=s<=N_E                             #
  #----------------------------------------------------------------------------------------------------------------------------#

  # Get all N's for which there is a max r (-N_C:(N_E-2)) for which P((E-C)<=r|Ha)<=beta+eps, and select that max r
  # Note that this does not take into account the first stage. However, the probability of rejecting Ha
  # must be lower in the second stage, compared to a 1-stage design, as some cases already rejected at first stage,
  # meaning that an rmax, calculated using a 1-stage design, is a good maximum
  #-----------------------------------------------------------------------------------------------------------------

  res0_t           <- data.frame(do.call("rbind", mapply(function(a,b,c) cbind(N=a,rtemp=(-c):(b-2),N_E=b,N_C=c),     # mapply outputs a list of matrices
                                                         a=SS,b=N_E,c=N_C,                     # do.call("rbind" ...) merges different
                                                         SIMPLIFY=F)))                                                       # list elements (one N-all r's) to one matrix


  res0_t$beta      <- mapply(function(a,b,c) bin_dif_cdf(z=a,n1=b,n2=c,p1=p_E,p2=p_C,type="normal"),
                                                         a=res0_t$rtemp,b=res0_t$N_E,c=res0_t$N_C,
                                                         SIMPLIFY=T)

  res0             <- aggregate(res0_t[res0_t$beta<=(beta +eps),]$rtemp, by = list(res0_t[res0_t$beta<=(beta +eps),c("N")  ],
                                                                                   res0_t[res0_t$beta<=(beta +eps),c("N_E")],
                                                                                   res0_t[res0_t$beta<=(beta +eps),c("N_C")]), FUN=max)
  names(res0 )     <- c("N","N_E","N_C","rmax")


  # Get for selected N all possible n1_C (1:(n_C-1)) and n1_E (n1_C*alloc)
  # + create r1max with -n1_C<=r1<=n1_E-1 and r1<=r-1 (note: r1<n1_E, otherwise first stage makes no sense: futility even if all outcomes a success)
  #------------------------------------------------------------------------------------------------------------------------------------------------

  res1           <- data.frame(do.call("rbind", mapply(function(a,b,c,d) cbind(N=a,N_E=b,N_C=c,n1_E=alloc*(1:(c-1)),n1_C=(1:(c-1)),rmax=d),
                                                       a=res0$N,b=res0$N_E,c=res0$N_C,d=res0$rmax,
                                                       SIMPLIFY=F)))
  res1$r1max     <- pmin(res1$rmax-1,res1$n1_E-1) # note: same difference in successes, but lower N--> cdf is higher, so rmax is a good maximum

  # Get for selected N,n1_E and n1_C --> r1's (-n1_C:[min(rmax,n1_E-1)], where P((x1_E-x1_C)<=r_1|Ha)<=beta+eps)
  # (note: r1<n1_E, otherwise first stage makes no sense: futility even if all outcomes a success)
  #-----------------------------------------------------------------------------------------------------

  res2_t      <- data.frame(do.call("rbind", mapply(function(a,b,c,d,e,f,g) cbind(N=a,N_E=b,N_C=c,n1_E=d,n1_C=e,r1=((-e):f),r1max=f,rmax=g),
                                                    a=res1$N,b=res1$N_E,c=res1$N_C,d=res1$n1_E,e=res1$n1_C,f=res1$r1max,g=res1$rmax,
                                                    SIMPLIFY=F)))
  res2_t$beta <- mapply(function(a,b,c) bin_dif_cdf(z=a,n1=b,n2=c,p1=p_E,p2=p_C,type="normal"),
                        a=res2_t$r1,b=res2_t$n1_E,c=res2_t$n1_C,
                        SIMPLIFY=T)
  res2        <- res2_t[res2_t$beta <= (beta +eps),]


  # Get for selected N, n1_E, n1_C and r1: r2's ((r1+1):rmax)
  #----------------------------------------------------------

  res3        <- data.frame(do.call("rbind", mapply(function(a,b,c,d,e,f,g) cbind(N=a,N_E=b,N_C=c,n1_E=d,n1_C=e,n2_E=b-d,n2_C=c-e,r1=g,r2=((g+1):f),rmax=f),
                                                    a=res2$N,b=res2$N_E,c=res2$N_C,d=res2$n1_E,e=res2$n1_C,f=res2$rmax,g=res2$r1,
                                                    SIMPLIFY=F)))

  #-------------------------------------------------------------------------#
  # Calculate beta and eta for all scenarios (only r1 and r2 needed; not s) #
  #-------------------------------------------------------------------------#

  # Split up in different datasteps, and 'break' loop  if beta condition not met
  res3_list        <- split(res3,list(res3$N,res3$n1_E,res3$n2_E,res3$r1),drop=T)
  pb <- tcltk::tkProgressBar(title = "Loop 1/4: Calculate beta", min = 0,max = length(res3_list), width = 300)

  for (i in 1:length(res3_list)){
    dt<-res3_list[[i]]
    dt$beta_temp<-NA
    tcltk::setTkProgressBar(pb, i)
    for (j in 1:dim(dt)[1]){
      dt[j,]$beta_temp <- probho_R(n1_E_=dt[j,]$n1_E,n1_C_=dt[j,]$n1_C,n2_E_=dt[j,]$n2_E,n2_C_=dt[j,]$n2_C,r1_=dt[j,]$r1,r_=dt[j,]$r2,p_E_=p_E,p_C_=p_C,type="normal")
      if (dt[j,]$beta_temp > (beta + eps) ) {
        break # need break to avoid unnecessary calculations
      }
    }
    if (i==1) {res3_list_new      <-list(dt)}   # Create list with first dataset
    if (i!=1) {res3_list_new[[i]] <-dt }        # Next iterations: append dataset
  }
  close(pb)
  res3             <- data.frame(do.call("rbind",res3_list_new))
  res3_1           <- res3[!is.na(res3$beta_temp) & res3$beta_temp <= (beta + eps),]

  # Split up in different datasteps, and 'break' loop  if eta condition not met
  res3_1_list        <- split(res3_1,list(res3_1$N,res3_1$n1_E,res3_1$n2_E,res3_1$r1),drop=T)
  pb <- tcltk::tkProgressBar(title = "Loop 2/4: Calculate eta", min = 0,max = length(res3_1_list), width = 300)

  for (i in 1:length(res3_1_list)){
    dt<-res3_1_list[[i]]
    dt$eta_temp<-NA
    tcltk::setTkProgressBar(pb, i)
    for (j in dim(dt)[1]:1){
      dt[j,]$eta_temp <- probho_R(n1_E_=dt[j,]$n1_E,n1_C_=dt[j,]$n1_C,n2_E_=dt[j,]$n2_E,n2_C_=dt[j,]$n2_C,r1_=dt[j,]$r1,r_=dt[j,]$r2,p_E_=p_M,p_C_=p_M,type="normal")
      if (dt[j,]$eta_temp < (eta - eps) ) {
        break # need break to avoid unnecessary calculations
      }
    }
    if (i==1) {res3_1_list_new      <-list(dt)}   # Create list with first dataset
    if (i!=1) {res3_1_list_new[[i]] <-dt }        # Next iterations: append dataset
  }
  close(pb)
  res3_1             <- data.frame(do.call("rbind",res3_1_list_new))
  res3_2             <- res3_1[!is.na(res3_1$eta_temp) & res3_1$eta_temp >= (eta - eps),]


  # Get for selected N,n1_E,n1_C,r1 and r2: s ((r2+2):N_E)
  #----------------------------------------------------------

  res4              <- data.frame(do.call("rbind", mapply(function(a,b,c,d,e,f,g,h,i,j,k) cbind(N=a,N_E=b,N_C=c,n1_E=d,n1_C=e,n2_E=f,n2_C=g,r1=h,r2=i,s=((i+2):b),beta_temp=j,eta_temp=k,alpha=alpha,beta=beta,eta=eta,pi=pi),
                                                          a=res3_2$N,b=res3_2$N_E,c=res3_2$N_C,d=res3_2$n1_E,e=res3_2$n1_C,f=res3_2$n2_E,g=res3_2$n2_C,h=res3_2$r1,i=res3_2$r2,j=res3_2$beta_temp,k=res3_2$eta_temp,
                                                          SIMPLIFY=F)))


  #----------------------------------------------------------#
  # Calculate pi and alpha for all scenarios (s needed)      #
  #----------------------------------------------------------#

  res4_list        <- split(res4,list(res4$N,res4$n1_E,res4$n2_E,res4$r1,res4$r2),drop=T)
  pb <- tcltk::tkProgressBar(title = "Loop 3/4: Calculate pi", min = 0,max = length(res4_list), width = 300)

  for (i in 1:length(res4_list)){
    dt<-res4_list[[i]]
    dt$pi_temp<-NA
    tcltk::setTkProgressBar(pb, i)
    for (j in 1:dim(dt)[1]){
      dt[j,]$pi_temp <- probha_R(n1_E_=dt[j,]$n1_E,n1_C_=dt[j,]$n1_C,n2_E_=dt[j,]$n2_E,n2_C_=dt[j,]$n2_C,r1_=dt[j,]$r1,s_=dt[j,]$s,p_E_=p_E,p_C_=p_C,type="normal")
      if (dt[j,]$pi_temp < (pi - eps)) {
        break # need break to avoid unnecessary calculations
      }
    }
    if (i==1) {res4_list_new      <-list(dt)}   # Create list with first dataset
    if (i!=1) {res4_list_new[[i]] <-dt }        # Next iterations: append dataset
  }
  close(pb)
  res4             <- data.frame(do.call("rbind",res4_list_new))
  res4_1           <- res4[!is.na(res4$pi_temp) & res4$pi_temp >= (pi - eps),]

  res4_1_list        <- split(res4_1,list(res4_1$N,res4_1$n1_E,res4_1$n2_E,res4_1$r1,res4_1$r2),drop=T)
  pb <- tcltk::tkProgressBar(title = "Loop 4/4: Calculate alpha", min = 0,max = length(res4_1_list), width = 300)

  for (i in 1:length(res4_1_list)){
    dt<-res4_1_list[[i]]
    dt$alpha_temp<-NA
    for (j in dim(dt)[1]:1){
      dt[j,]$alpha_temp <- probha_R(n1_E_=dt[j,]$n1_E,n1_C_=dt[j,]$n1_C,n2_E_=dt[j,]$n2_E,n2_C_=dt[j,]$n2_C,r1_=dt[j,]$r1,s_=dt[j,]$s,p_E_=p_M,p_C_=p_M,...)
      if (dt[j,]$alpha_temp > (alpha + eps)) {
        break # need break to avoid unnecessary calculations
      }
    }
    if (i==1) {res4_1_list_new      <-list(dt)}   # Create list with first dataset
    if (i!=1) {res4_1_list_new[[i]] <-dt }        # Next iterations: append dataset
  }
  close(pb)
  res4_1           <- data.frame(do.call("rbind",res4_1_list_new))
  res5             <- res4_1[!is.na(res4_1$alpha_temp) & res4_1$alpha_temp <= (alpha + eps),]

  res5$PET.p0      <- mapply(function (a,b,c,d,e,f) bin_dif_cdf(z=c,n1=a,n2=b,p1=p_C,p2=p_C,...),
                             a=res5$n1_E,b=res5$n1_C,c=res5$r1,
                             SIMPLIFY=T)
  res5$EN.p0       <- res5$N-(res5$n2_E+res5$n2_C)*res5$PET.p0

  res5             <- data.table::as.data.table(res5)

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

  res$lambda<-1-(res$eta +res$alpha)
  res$delta <-1-(res$beta+res$pi)

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

  res$n1<-res$n1_E+res$n1_C
  res$n2<-res$n2_E+res$n2_C
  res$alloc<-alloc
  res<- res[,-c("N_min","EN.p0_min","EN.p0_N_min")]

  if (alloc==1) {
    res$d0<-round(100*(res$r1/res$n1_E))
    res$d1<-round(100*(res$r2/res$N_E))
    res$d2<-round(100*(res$s /res$N_E))
    res<- cbind(res[,c("r1","d0","n1_E","n1_C","n1","r2","d0","s","d2","n2_E","n2_C","n2","N","N_E","N_C","EN.p0","PET.p0","MIN","OPT","ADMISS","alpha","beta","eta","pi","lambda","delta","alloc")],p0=p0,pa=pa,res[,c("alpha_param","beta_param","eta_param","pi_param")])
  }

  if (alloc>1) {
    res<- cbind(res[,c("r1",     "n1_E","n1_C","n1","r2"     ,"s"     ,"n2_E","n2_C","n2","N","N_E","N_C","EN.p0","PET.p0","MIN","OPT","ADMISS","alpha","beta","eta","pi","lambda","delta","alloc")],p0=p0,pa=pa,res[,c("alpha_param","beta_param","eta_param","pi_param")])
  }

  res<-res[order(res$N,res$n1_E,res$n2_E,res$r1,res$r2,res$s),]
  return(res)
}


# TEST (Sargent DJ, Goldberg RM. A Three-Outcome Design for Phase II Clinical Trials. Controlled Clinical Trials 22:117-125
#--------------------------------------------------------------------------------------------------------------------------
# test_0_a<- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,2),pa=rep((a+0.2 ),2),alpha=c(0.05,0.1),beta=c(0.1,0.1),eta=c(0.85,0.8),pi=c(0.8,0.8)),a=c(0.1,0.2,0.3,0.4),SIMPLIFY=F)))
# test_0_b<- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,2),pa=rep((a+0.15),2),alpha=c(0.05,0.1),beta=c(0.1,0.1),eta=c(0.85,0.8),pi=c(0.8,0.8)),a=c(0.1,0.2,0.3,0.4),SIMPLIFY=F)))
#
# test_0<-cbind(rbind(test_0_a,test_0_b),data.frame(N_min=c(38+10,28+10,52+11,40+13,59+10,43+13,59+10,47+14,65+10,46+14,75+31,72+10,64+64,30+67,91+46,91+10)*2-10,
#                                                   N_max=c(19+32,20+20,32+36,24+36,30+44,36+23,36+47,24+41,33+45,26+37,50+65,35+53,52+78,45+58,52+89,55+56)*2+10))
#
# for (i in 1:dim(test_0)[1]){
#   print(paste0("test=",i,"/",dim(test_0)[1]))
#   res<- sargent2stage_R(p0=test_0[i,]$p0,pa=test_0[i,]$pa,alpha=test_0[i,]$alpha,beta=test_0[i,]$beta,eta=test_0[i,]$eta,pi=test_0[i,]$pi,
#                                     eps=0.005,N_min=test_0[i,]$N_min,N_max=test_0[i,]$N_max,alloc=1,type="normal")
#
#   if (i==1) {test_list_O      <-list(res[OPT=="Optimal",])}   # Create list with first dataset
#   if (i!=1) {test_list_O[[i]] <-res[OPT=="Optimal",] }        # Next iterations: append dataset
#
#   if (i==1) {test_list_M      <-list(res[MIN=="Minimax",])}   # Create list with first dataset
#   if (i!=1) {test_list_M[[i]] <-res[MIN=="Minimax",] }        # Next iterations: append dataset
#
# }
# test_O<- data.frame(do.call("rbind",test_list_O))
# test_O<-test_O[,c("p0","pa","alpha_param","beta_param","eta_param","pi_param","alpha","beta","eta","pi","lambda","delta","n1_E","n1_C","n2_E","n2_C","n1","n2","r1","r2","s","N","EN.p0")]
#
# test_M<- data.frame(do.call("rbind",test_list_M))
# test_M<-test_M[,c("p0","pa","alpha_param","beta_param","eta_param","pi_param","alpha","beta","eta","pi","lambda","delta","n1_E","n1_C","n2_E","n2_C","n1","n2","r1","r2","s","N","EN.p0")]
