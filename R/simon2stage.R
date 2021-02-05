
P_Simon_reject_Ha <- function(n1, n2, r, b_p0, B_p0, b_pa, B_pa){
  ## r is the total number of seen cases on both n1 and n2 together
  ## Calculate alpha/beta for all values up to r (see formula (1) paper Simon)

  r_1 <- (max(0,r-n2)):(max(0,min(n1-1, r-1))) # for min value: r_1 must be minimal (r-n2): if n2<r, then x1>=(r+1-n2), so r1>=r-n2
                                               # for max value: r_1<n1, otherwise first stage makes no sense: futility even if all outcomes a success
  x1  <- r_1 + 1  # e.g. when r=2; r_1=0 or r_1=1; x1=1 or x1=2

  # alpha
  B_p0_r1    <- B_p0[[n1]][1 + r_1]     # P(X1<=r),    so pbinom(r_1 ,n1,p0), indexing with "+1", as B_p0 starts from 0 successes
  b_p0_x1    <- b_p0[[n1]][1 + x1]      # P(X1=r) ,    so dbinom(x1  ,n1,p0), indexing with "+1", as b_p0 starts from 0 successes
  B_p0_r2    <- B_p0[[n2]][1 + (r-x1)]  # P(X1+X2<=r), so pbinom(r-x1,n2,p0), indexing with "+1", as B_p0 starts from 0 successes
  eta_temp   <- B_p0_r1 + rev(cumsum(rev(b_p0_x1 * B_p0_r2))) # e.g. for r=2: r_1=0: P(X1<=0) + P(X1=1)P(X1+X2<=2) + P(X1=2)P(X1+X2<=2)
                                                             #                r_1=1: P(X1<=1) +                      P(X1=2)P(X1+X2<=2)
  # beta
  B_pa_r1    <- B_pa[[n1]][1 + r_1]
  b_pa_x1    <- b_pa[[n1]][1 + x1]
  B_pa_r2    <- B_pa[[n2]][1 + (r-x1)]
  beta_temp  <- B_pa_r1 + rev(cumsum(rev(b_pa_x1 * B_pa_r2)))
  list(N = n1 + n2, n1 = n1, r1 = r_1, n2 = n2, r2 = r, alpha_temp = 1 - eta_temp, beta_temp = beta_temp)
}

#' @title Simon 2-stage function
#' @description
#' The primary objective of a phase II clinical trial of a new drug or regimen is to determine whether it has sufficient biological activity
#' against the disease under study to warrant more extensive development. \cr
#' This function calculates the sample size needed in a Simon 2-stage design which is a
#' two-stage design that is optimal in the sense that the expected sample size is minimized if the
#' regimen has low activity subject to constraints upon the size of the type 1 and type 2 errors. \cr
#' Two-stage designs which minimize the maximum sample size are also determined.
#' @param p0 probability of the uninteresting response (null hypothesis \eqn{H0})
#' @param pa probability of the interesting response (alternative hypothesis Ha)
#' @param alpha Type I error rate \eqn{P(reject H0|H0)}
#' @param beta Type II error rate \eqn{P(reject Ha|Ha)}
#' @param eps tolerance default value = 0.005
#' @param N_min minimum sample size value for grid search
#' @param N_max maximum sample size value for grid search
#' @param int pre-specified interim analysis percentage information
#' @param int_window window around interim analysis percentage (e.g. 0.5 +- 0.025). 0.025 is default value
#' @return a data.frame with elements
#' \itemize{
#' \item n1: total number of patients in stage1
#' \item n2: total number of patients in stage2
#' \item N: total number of patients=n1+n2
#' \item r1: ("r" stands for "rejection") threshold for "rejecting" Ha: if x1<=r1 --> stop for futility at first stage
#' \item r2: ("r" stands for "rejection") threshold for "rejecting" Ha: if x1+x2<=r --> futility at second stage
#' \item eff: (r2 + 1)/N
#' \item 90%CI_low: Result of call to OneArmPhaseTwoStudy::get_CI. Confidence interval according to Koyama and Chen (1989)
#' \item 90%CI_high: Result of call to OneArmPhaseTwoStudy::get_CI. Confidence interval according to Koyama and Chen (1989)
#' \item EN.p0: expected sample size under H0
#' \item PET.p0: probability of terminating the trial at the end of the first stage under H0
#' \item MIN: column indicating if the design is the minimal design
#' \item OPT: column indicating if the setting is the optimal design
#' \item ADMISS: column indicating if the setting is the admissible design
#' \item alpha: the actual alpha value which is smaller than \code{alpha_param + eps}
#' \item beta: the actual beta value where which is smaller than \code{beta_param + eps}
#' \item p0: your provided \code{p0} value
#' \item pa: your provided \code{pa} value
#' \item alpha_param: your provided \code{alpha} value
#' \item beta_param: your provided \code{beta} value
#' }
#' @details
#' if x1<=r1 --> stop futility \cr
#' if (x1+x2)<=r --> futility \cr
#' if (x1+x2)> r --> efficacy \cr
#' @references Simon R. Optimal two-stage designs for phase II clinical trials. Control Clin Trials. 1989;10(1):1-10. doi:10.1016/0197-2456(89)90015-9
#'     Koyama T, Chen H. Proper inference from simon’s two-stage designs. Stat Med. 2008; 27:3145–154;
#' @export
#' @examples
#' samplesize <- simon2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.2,
#'                           eps = 0.005, N_min = 1, N_max = 50)
#' plot(samplesize)
#' \donttest{
#' samplesize <- simon2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.2,
#'                           eps = 0.005, N_min = 1, N_max = 50,int=0.33,int_window=0.025)
#' plot(samplesize)
#' }
#' \donttest{
#' ## Example 1
#' test <- data.frame(p0 = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
#'                    pa = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7) + 0.2)
#' test <- merge(test,
#'               data.frame(alpha = c(0.1, 0.05, 0.05), beta = c(0.1, 0.2, 0.1)))
#' samplesize <- fleming1stage(p0 = test$p0, pa = test$pa, alpha = test$alpha, beta = test$beta)
#' samplesize <- simon2stage(p0 = test$p0, pa = test$pa, alpha = test$alpha, beta = test$beta,
#'                           N_min = 10, N_max = samplesize$n + 15)
#' optimal_minimax <- lapply(samplesize, FUN=function(x){
#'   cbind(subset(x, OPT == "Optimal", c("r1","n1","r2","N","EN.p0","PET.p0")),
#'         subset(x, MIN == "Minimax", c("r1","n1","r2","N","EN.p0","PET.p0")))
#' })
#' optimal_minimax <- data.table::rbindlist(optimal_minimax)
#'
#' ## Example 2
#' test <- data.frame(p0 = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
#'                    pa = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7) + 0.15)
#' test <- merge(test,
#'               data.frame(alpha = c(0.1, 0.05, 0.05), beta = c(0.1, 0.2, 0.1)))
#' samplesize <- fleming1stage(p0 = test$p0, pa = test$pa, alpha = test$alpha, beta = test$beta)
#' samplesize <- simon2stage(p0 = test$p0, pa = test$pa, alpha = test$alpha, beta = test$beta,
#'                           N_min = 25, N_max = samplesize$n + 20)
#' optimal_minimax <- lapply(samplesize, FUN=function(x){
#'   cbind(subset(x, OPT == "Optimal", c("r1","n1","r2","N","EN.p0","PET.p0")),
#'         subset(x, MIN == "Minimax", c("r1","n1","r2","N","EN.p0","PET.p0")))
#' })
#' optimal_minimax <- data.table::rbindlist(optimal_minimax)
#' }

simon2stage <- function(p0, pa, alpha, beta, eps = 0, N_min, N_max, int=0, int_window=0.025){

  if(length(p0) > 1 && length(pa) > 1){
    results <- mapply(null = p0, alternative = pa, alpha = alpha, beta = beta, eps = eps, N_min = N_min, N_max = N_max, int = int, int_window = int_window,
                      FUN = function(null, alternative, alpha, beta, eps, N_min, N_max, int, int_window){
                        simon2stage.default(p0 = null, pa = alternative, alpha = alpha, beta = beta, eps = eps, N_min = N_min, N_max = N_max, int = int, int_window = int_window)
                      }, SIMPLIFY = FALSE)
  }else{
    results <- simon2stage.default(p0 = p0, pa = pa, alpha = alpha, beta = beta, eps = eps, N_min = N_min, N_max = N_max, int = int, int_window = int_window)
  }
}

simon2stage.default <- function(p0, pa, alpha, beta, eps = 0, N_min, N_max, int=0, int_window=0.025) {

  # WARNING MESSAGES
  if (pa < p0) {
    stop("p0 should be smaller than pa")
  }

  if (N_min <1 ) {
    stop("N_min should be >=1")
  }

  EN.p0 <- EN.p0_N_min <- EN.p0_N_n1_min <- EN.p0_min <- EN.p0_N_min_int <- MIN <- N <- OPT <- NULL
  rowid <- n1 <- n2 <- r <- r2 <- NULL

  #----------------------------------------------------------------------------------------------------------#
  # Get all possible scenarios for N, n1, r1, r2 and n2                                                      #
  # Note that this is the possible range of values: 0 <=r1<n1; r1+1<=r; r<=n1+n2                             #
  #----------------------------------------------------------------------------------------------------------#

  # Get all N's for which there is a max r (1:N) for which P(X<=r|Ha)<=beta+eps, and select that max r
  # Note that this is not taking into account first stage. However, the probability of rejecting Ha
  # should be lower in the second stage, compared to a 1-stage design, as some cases already rejected at first stage,
  # meaning that an rmax, calculated using a 1-stage design, is a good maximum
  #------------------------------------------------------------------------------------------------------------

  res0 <- lapply(N_min:N_max, FUN=function(N) cbind(N = N, rtemp = 1:N))
  res0 <- do.call(rbind, res0)
  res0 <- data.frame(res0)
  res0$betamax <- pbinom(q = res0$rtemp, size = res0$N, prob = pa, lower.tail = TRUE)
  res0 <- res0[!is.na(res0$betamax) & res0$betamax <= (beta + eps), ]
  # Select all possible N"s: there needs to be
  # at least one r with P(X<=r|Ha)<=beta+eps
  res0 <- aggregate(res0$rtemp, by = list(res0$N), max)
  names(res0) <- c("N", "rmax")

  # Get for selected N's all possible n1's (1:N-1)
  #-----------------------------------------------
  res1 <- mapply(N = res0$N, rmax = res0$rmax, FUN = function(N, rmax) cbind(N = N, n1 = (1:(N - 1)), rmax = rmax), SIMPLIFY = FALSE)
  res1 <- do.call(rbind, res1)
  res1 <- data.frame(res1)

  # Get for selected N's and n1's: r's (1:rmax)
  #--------------------------------------------
  res2 <- mapply(N = res1$N, n1 = res1$n1, rmax = res1$rmax,
                 FUN = function(N, n1, rmax) cbind(N = N, n1 = n1, rmax = rmax, r = (1:rmax)),
                 SIMPLIFY = FALSE)
  res2 <- do.call(rbind, res2)
  res2 <- data.frame(res2)
  res2$n2 <- res2$N - res2$n1

  # Calculate beta and alpha
  #-------------------------

  ###   use a lookup table for formula (1) from paper Simon
  ###   b: probability mass       (see formula (1) paper Simon)
  ###   B: cumulative probability (see formula (1) paper Simon)
  nmax <- N_max
  b_p0 <- lapply(1:nmax, FUN = function(n) dbinom(0:nmax, size = n, prob = p0))
  b_pa <- lapply(1:nmax, FUN = function(n) dbinom(0:nmax, size = n, prob = pa))
  B_p0 <- lapply(1:nmax, FUN = function(n) pbinom(0:nmax, size = n, prob = p0, lower.tail = TRUE))
  B_pa <- lapply(1:nmax, FUN = function(n) pbinom(0:nmax, size = n, prob = pa, lower.tail = TRUE))

  res2 <- data.table::setDT(res2)
  res2$rowid <- seq_len(nrow(res2))
  res3 <- res2[, P_Simon_reject_Ha(n1 = n1, n2 = n2, r = r, b_p0 = b_p0, B_p0 = B_p0, b_pa = b_pa, B_pa = B_pa), by = list(rowid)]
  res3 <- data.table::setDF(res3)
  res3$diff_beta <- res3$beta_temp - beta
  res4 <- res3[res3$diff_beta <= eps, ]
  res4$diff_alpha <- res4$alpha_temp - alpha
  res5 <- res4[res4$diff_alpha <= eps, ]

  res5[2:6] <- lapply(res5[2:6], as.numeric) # convert integer into numeric

  if(nrow(res5) == 0){
    stop("No data satisfying the H0/Ha criteria")
  }

  res5$alpha <- alpha
  res5$beta <- beta

  res5$PET.p0 <- pbinom(q = res5$r1, size = res5$n1, prob = p0, lower.tail = TRUE)
  res5$EN.p0 <- res5$N - ((res5$N - res5$n1) * res5$PET.p0)

  res5 <- data.table::setDT(res5)

  res5 <- res5[, N_min          := min(N)]
  res5 <- res5[, EN.p0_min      := min(EN.p0)]
  res5 <- res5[, EN.p0_N_min    := min(EN.p0) , by = list(N)]
  res5 <- res5[, EN.p0_N_n1_min := min(EN.p0) , by = list(N,n1)]

  if (int>0){
    res6_int <- res5[EN.p0 == EN.p0_N_n1_min & floor((int-int_window)*N)<=n1 & n1<=ceiling((int+int_window)*N), ]
    res6_int <- res6_int[, EN.p0_N_min_int := min(EN.p0) , by = list(N)]
    res6_int <- res6_int[EN.p0 == EN.p0_N_min_int   , ]
    res6_int$INTERIM=int;
  }

  res6         <- res5[EN.p0 == EN.p0_N_min   , ]
  res6$OPT     <- res6$MIN <- res6$ADMISS <- c("")
  res6$OPT[which(res6$EN.p0 == res6$EN.p0_min)] <- "Optimal"
  res6$MIN[which(res6$N == res6$N_min)] <- "Minimax" # Note: if multiple designs that meet the criteria for minimax:choose
  # design with minimal expected sample size under H0: "Optimal minimax design"

  res6 <- data.table::setDF(res6)

  # Get admissible designs
  y <- res6[, c("N", "EN.p0")]
  con.ind <- chull(y)[chull((y)) == cummin(chull((y)))]
  res <- res6[res6$N >= min(res6$N[res6$MIN == "Minimax"]) & res6$N <= max(res6$N[res6$OPT == "Optimal"]), ]
  res$ADMISS[which((rownames(res) %in% c(con.ind)) & (res$N > res$N_min) & (res$EN.p0 > res$EN.p0_min) & (res$N < res$N[res$OPT == "Optimal"]))] <- "Admissible"

  if (int>0){
    res$EN_opt  <- "EN.p0_optimal"
    res6_int<-res6_int[res6_int$N >= min(res6$N[res6$MIN == "Minimax"]) & res6_int$N <= max(res6$N[res6$OPT == "Optimal"]), ]
    res<-merge(x = res, y = res6_int, by = c("N","n1","n2","r1","r2","alpha_temp","beta_temp","diff_beta","diff_alpha","alpha","beta","PET.p0","EN.p0","N_min","EN.p0_min","EN.p0_N_min","EN.p0_N_n1_min","rowid"), all = TRUE)

    res$ADMISS[is.na(res$ADMISS)] <- "" # AFter merging: replace NA's with """
    res$MIN   [is.na(res$MIN)   ] <- ""
    res$OPT   [is.na(res$OPT)   ] <- ""
    res$EN_opt[is.na(res$EN_opt)] <- ""
  }


  # Calculate 1-2*alpha confidence interval, based on Koyama, Statistics in Medicine 2008

  res$eff <- paste0(res$r2 + 1, "/", res$N, " (", 100 * round((res$r2 + 1) / res$N, 3), "%)")
  CI <- mapply(a = res$r2 + 1, b = res$r1, c = res$n1, d = res$N,
               FUN = function(a, b, c, d) OneArmPhaseTwoStudy::get_CI(k = a, r1 = b, n1 = c, n = d, alpha = alpha, precision = 3))
  res$CI_low  <- 100 * unlist(CI[rownames(CI) == "CI_low", ])
  res$CI_high <- 100 * unlist(CI[rownames(CI) == "CI_high", ])

  res <- data.table::setnames(res,
                              old = c("alpha", "beta"),
                              new = c("alpha_param", "beta_param"))
  res <- data.table::setnames(res,
                              old = c("alpha_temp", "beta_temp"),
                              new = c("alpha", "beta"))

  if (int>0){
  res <- cbind(design_nr=1:dim(res)[1],
               res[, c("r1", "n1", "r2", "n2", "N", "eff", "CI_low", "CI_high", "EN.p0", "PET.p0", "MIN", "OPT", "ADMISS","EN_opt","INTERIM", "alpha", "beta")],
               p0 = p0, pa = pa,
               res[, c("alpha_param", "beta_param")])
  }
  if (int==0){
    res <- cbind(design_nr=1:dim(res)[1],
                 res[, c("r1", "n1", "r2", "n2", "N", "eff", "CI_low", "CI_high", "EN.p0", "PET.p0", "MIN", "OPT", "ADMISS","alpha", "beta")],
                 p0 = p0, pa = pa,
                 res[, c("alpha_param", "beta_param")])
  }
  res <- data.table::setnames(res,
                              old = c("CI_low", "CI_high"),
                              new = c(paste0(100 - 2 * 100 * alpha, "%CI_low"), paste0(100 - 2 * 100 * alpha, "%CI_high")))
  attr(res, "inputs") <- list(p0 = p0, pa = pa, alpha = alpha, beta = beta, eps = eps, N_min = N_min, N_max = N_max)
  class(res) <- c("2stage", "simon", "data.frame")
  res
}



# TEST (Simon R. Optimal Two-Stage Designs for Phase II Clinical Trials. Controlled Clinical Trials 1989;10:1-10 : Table 1 and 2
#-------------------------------------------------------------------------------------------------------------------------------
# test_1_0  <- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,3),pa=rep((a+0.2),3),alpha=c(0.1,0.05,0.05),beta=c(0.1,0.2,0.1)),a=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7),SIMPLIFY=F)))
# for(i in 1:dim(test_1_0)[1]){
#  nmax   <- PhIIdesign::fleming1stage(p0=test_1_0[i,]$p0,pa=test_1_0[i,]$pa,alpha=test_1_0[i,]$alpha,beta=test_1_0[i,]$beta,eps=0)$n+15
#  res    <- simon2stage  (p0=test_1_0[i,]$p0,pa=test_1_0[i,]$pa,alpha=test_1_0[i,]$alpha,beta=test_1_0[i,]$beta,eps=0,N_min=10,N_max=nmax)
#  res_dt <- cbind(res[OPT=="Optimal",c("r1","n1","r2","N","EN.p0","PET.p0")],res[MIN=="Minimax",c("r1","n1","r2","N","EN.p0","PET.p0")])
#  if (i==1) {test1_list      <- list(res_dt)}
#  if (i!=1) {test1_list[[i]] <- res_dt }
# }
# test1             <- cbind(test_1_0,data.frame(do.call("rbind",test1_list)))
#
# test_2_0  <- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,3),pa=rep((a+0.15),3),alpha=c(0.1,0.05,0.05),beta=c(0.1,0.2,0.1)),a=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),SIMPLIFY=F)))
# for(i in 1:dim(test_2_0)[1]){
#  nmax   <- PhIIdesign::fleming1stage(p0=test_2_0[i,]$p0,pa=test_2_0[i,]$pa,alpha=test_2_0[i,]$alpha,beta=test_2_0[i,]$beta,eps=0)$n+20
#  res    <- PhIIdesign::simon2stage  (p0=test_2_0[i,]$p0,pa=test_2_0[i,]$pa,alpha=test_2_0[i,]$alpha,beta=test_2_0[i,]$beta,eps=0,N_min=25,N_max=nmax)
#  res_dt <- cbind(res[OPT=="Optimal",c("r1","n1","r2","N","EN.p0","PET.p0")],res[MIN=="Minimax",c("r1","n1","r2","N","EN.p0","PET.p0")])
#  if (i==1) {test2_list      <- list(res_dt)}
#  if (i!=1) {test2_list[[i]] <- res_dt }
# }
# test2             <- cbind(test_2_0,data.frame(do.call("rbind",test2_list)))
