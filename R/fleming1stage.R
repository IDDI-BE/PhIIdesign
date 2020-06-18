
#' @title The Fleming 1-stage function
#' @description This function calculates sample sizes of the Fleming 1-stage design, for p0<pa
#' @details n= total number of patients
#' @details if x<=r --> futility
#'
#' @param p0 uninteresting response (null hypothesis H0)
#' @param pa interesting response (alternative hypothesis Ha)
#' @param alpha P(reject H0|H0)
#' @param beta P(reject Ha|Ha)
#' @param eps tolerance (actual alpha<=alpha+eps; actual beta<=beta+eps; actual eta>=eta-eps; actual p>=pi-eps); default value = 0.005
#' @examples
#' fleming1stage(p0=0.45,pa=0.7,alpha=0.05,beta=0.2)
#' @export
#' @importFrom stats pbinom qbinom

#p0=0.1;pa=0.3;alpha=0.05;beta=0.1;eps=0.005

fleming1stage<-function (p0,pa,alpha,beta,eps=0.005)
{
  if (pa<p0) {stop('p0 should be smaller than pa')}

  n         <- 0
  beta_temp <- 1

  while (beta_temp > beta+eps) {
    n         <- n + 1
    r_temp    <- qbinom(p=1-(alpha+eps),size=n,prob=p0,lower.tail=T)
    beta_temp <- pbinom(q=r_temp       ,size=n,prob=pa,lower.tail=T)
  }

  alpha_temp= 1-pbinom(q=r_temp,size=n,prob=p0,lower.tail=T)
  res<-data.frame(n=n,r=r_temp,alpha=alpha_temp,beta=beta_temp,p0=p0,pa=pa,alpha_param=alpha,beta_param=beta)
  return(res)
}

# TEST (A'Hern RP. Sample size tables for exact single-stage phase II designs. Statistics in Medicine 2001;20:859-866: Table 1
#------------------------------------------------------------------------------------------------------------------------------
# test0                    <- data.frame(do.call("rbind", mapply(function(a) cbind(p0=a,pa=seq(a,0.95,by=0.05)),a=seq(0.05,0.95,by=0.05),SIMPLIFY=F)))
# test                     <- test0[test0$pa>test0$p0,]
# test_alpha_0.05_beta_0.2 <- cbind(test,data.frame(do.call("rbind",mapply(function (a,b) PhIIdesign::fleming1stage(p0=a,pa=b,alpha=0.05,beta=0.2,eps=0.0005),a=test$p0,b=test$pa,SIMPLIFY=F))))
# test_alpha_0.05_beta_0.1 <- cbind(test,data.frame(do.call("rbind",mapply(function (a,b) PhIIdesign::fleming1stage(p0=a,pa=b,alpha=0.05,beta=0.1,eps=0.0005),a=test$p0,b=test$pa,SIMPLIFY=F))))
# test_alpha_0.01_beta_0.2 <- cbind(test,data.frame(do.call("rbind",mapply(function (a,b) PhIIdesign::fleming1stage(p0=a,pa=b,alpha=0.01,beta=0.2,eps=0.0005),a=test$p0,b=test$pa,SIMPLIFY=F))))
# test_alpha_0.01_beta_0.1 <- cbind(test,data.frame(do.call("rbind",mapply(function (a,b) PhIIdesign::fleming1stage(p0=a,pa=b,alpha=0.01,beta=0.1,eps=0.0005),a=test$p0,b=test$pa,SIMPLIFY=F))))
