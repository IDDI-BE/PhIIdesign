
#' @title Onestage combination function for optimal designs
#' @description This function calculates sample sizes of the Fleming 1-stage , Simon-2stage, Sargent's 1- and 2-stage design. calculations are performed jointly for different sets of proportions (p0, pA) and operating characteristics. All designs are summarized in 2 tables, 1 for one-stage designs (Fleming and sargent) and 1 for two-stage designs (Simon and sargent)
#' @param p0 uninteresting response (null hypothesis H0), can be a vector
#' @param pa interesting response (alternative hypothesis Ha), can be a vector, always same length as p0. The corresponding elements of p0 and pa are taken as a set of proportions
#' @param alpha P(reject H0|H0) for Fleming and Simon designs, can be a vector
#' @param beta P(reject Ha|Ha)  for Fleming and Simon designs, can be a vector
#' @param alpha2 P(reject H0|H0) for Sargent designs, can be a vector, same lenth as alpha
#' @param beta2 P(reject Ha|Ha)  for Sargent designs, can be a vector, same length as beta
#' @param eta P(reject Ha|H0) for Sargent designs, can be a vector, same length as alpha2
#' @param pi P(reject H0|Ha) for Sargent designs, can be a vector, same length as beta2
#' @param eps tolerance (actual alpha<=alpha+eps; actual beta<=beta+eps; actual eta>=eta-eps; actual p>=pi-eps); default value = 0.005
#' @export
#' @examples
#' allsinglearm(p0 = 0.1, pa = 0.7,
#'              alpha = 0.05, beta = 0.2, beta2 = 0.1, pi = 0.8, eta = 0.8)
allsinglearm<-function (p0,pa,alpha,beta,alpha2=alpha,beta2,pi,eta,eps=0.005){

onestage<-NULL
twostage<-NULL

for (p in 1:length(p0)){

  for (a in 1:length(alpha)){

    for (b in 1:length(beta)){

      #Get designs
      #-----------------------
      flem<-fleming1stage(p0=p0[p],pa=pa[p], alpha=alpha[a], beta=beta[b], eps = eps)
      sim<-simon2stage(p0=p0[p],pa=pa[p], alpha=alpha[a], beta=beta[b], eps = eps,
                      N_min=max(1,(flem$N-5)), N_max=(flem$N+15))
      sim<-(subset(sim,sim$OPT=="Optimal"))
      sar1<-(sargent1stage(p0=p0[p],pa=pa[p], alpha=alpha2[a], beta=beta2[b],
                          pi=pi[b],eta=eta[a],eps = eps, N_min=max(3,(flem$N-15)), N_max=(flem$N+15)))
      sar2<-sargent2stage(p0=p0[p],pa=pa[p], alpha=alpha2[a], beta=beta2[b],
                         pi=pi[b],eta=eta[a],eps = eps, N_min=max(3,(sim$N-15)), N_max=(sim$N+15))
      sar2<-subset(sar2,sar2$OPT=="Optimal")

      #Create overview tables
      #---------------------

      one<-c(p0[p],pa[p],alpha[a],beta[b],unlist(flem[1,1:2]),
            alpha2[a],beta2[b],pi[b],eta[a],unlist(sar1[1,3:1]))
      onestage<-rbind(onestage,one)

      two<-c(p0[p],pa[p],alpha[a],beta[b],unlist(sim[1,1:5]),
            alpha2[a],beta2[b],pi[b],eta[a],unlist(sar2[1,1:6]))
      twostage<-rbind(twostage,two)

      rm(one, two,flem,sim,sar1,sar2)

    }
  }
}

colnames(onestage)<-c("p0","pa","alpha","beta","N_Flem","R_Flem",
                   "alpha","beta","pi","eta","N_Sar","S_Sar","R_Sar")
colnames(twostage)<-c("p0","pa","alpha","beta","R1_Sim","N1_Sim","R_sim","N2_Sim","N_Sim",
                   "alpha","beta","pi","eta","R1_Sar","N1_Sar","R_Sar","S_Sar","N2_Sar","N_Sar")
res<-list(onestage,twostage)
return(res)

}
