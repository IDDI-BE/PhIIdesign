
#-----------------------------------------------------#
# two error probabilities functions (Equations p3530) #
#-----------------------------------------------------#

P_Sargent_reject_Ha <- function(n1, n2, r, b_p0, B_p0, b_pa, B_pa){
  ## r is the total number of seen cases on both n1 and n2 together
  ## Calculate eta/beta for all values up to r (see formula (1) paper Simon)

  r_1 <- (max(0,r-n2)):(max(0,min(n1-1, r-1))) # for min value: r_1 must be minimal (r-n2): if n2<r
  #(note: r_1<n1, otherwise first stage makes no sense: futility even if all outcomes a success)
  x1  <- r_1 + 1  # e.g. when r=2; r_1=0 or r_1=1; x1=1 or x1=2

  # eta
  B_p0_r1    <- B_p0[[n1]][1 + r_1]    # P(X1<=r_1),  so pbinom(r_1 ,n1,p0), indexing with "+1", as B_p0 starts from 0 successes
  b_p0_x1    <- b_p0[[n1]][1 + x1]     # P(X1=x1) ,   so dbinom(x1  ,n1,p0), indexing with "+1", as b_p0 starts from 0 successes
  B_p0_r2    <- B_p0[[n2]][1 + (r-x1)] # P(X1+X2>=r), so pbinom(r-x1,n2,p0), indexing with "+1", as B_p0 starts from 0 successes
  eta_temp   <- B_p0_r1 + rev(cumsum(rev(b_p0_x1 * B_p0_r2))) # e.g. for r=2: r_1=0: P(X1<=0) + P(X1=1)P(X1+X2<=2) + P(X1=2)P(X1+X2<=2)
                                                              #               r_1=1: P(X1<=1) +                      P(X1=2)P(X1+X2<=2)
  # beta
  B_pa_r1    <- B_pa[[n1]][1 + r_1]
  b_pa_x1    <- b_pa[[n1]][1 + x1]
  B_pa_r2    <- B_pa[[n2]][1 + (r-x1)]
  beta_temp  <- B_pa_r1 + rev(cumsum(rev(b_pa_x1 * B_pa_r2)))
  list(N = n1 + n2, n1 = n1, r1 = r_1, n2 = n2, r2 = r, eta_temp = eta_temp, beta_temp = beta_temp)
}

P_Sargent_reject_H0 <- function(n1, n2, r1min, s, b_p0, B_p0, b_pa, B_pa){
  ## Equation 4 page 122 Sargent et. Al

  r_1 <- r1min:(min(n1-1, s-1))
  x1  <- r_1 + 1

  # alpha
  B_p0_r1    <- (1-B_p0[[n1]][1 + s]) * as.numeric(n1 > s)  # P(X1>=s+1) =1-P(X1<=s)                   , so 1-pbinom(s     ,n1,p0,lower.tail=F), indexing with "+1", as B_p0_ut starts from 0 successes
  b_p0_x1    <-    b_p0[[n1]][1 + x1]                       # P(X1=x1)                                 , so   dbinom(x1    ,n1,p0)
  B_p0_s     <-  1-B_p0[[n2]][1 + (s-x1)]                   # P(X1+X2>=s+1)=1-P(X1+X2<=s)=1-P(X2<=s-x1), so 1-pbinom(s-x1-1,n2,p0)
  b_p0_s     <-    b_p0[[n2]][1 + (s-x1)]                   # P(X2=s)                                  , so   dbinom(s     ,n2,p0)
                                                            # B_p0_s + b_p0_s = P(X1+X2>=s). Note that calculating P(X1+X2>=s)=1-P(X1+X2<=s-1)
                                                            # is not possible here: if x1=s, then in P(X2<=s-x1-1), right side is negative, so no index
                                                            # in listings with binomial results
  alpha_temp <- B_p0_r1 + rev(cumsum(rev(b_p0_x1 * (b_p0_s+B_p0_s))))

  # pi
  B_pa_r1    <- (1-B_pa[[n1]][1 + s]) * as.numeric(n1 > s)  # P(X1>=s+1) =1-P(X1<=s)                   , so 1-pbinom(s     ,n1,pa,lower.tail=F), indexing with "+1", as B_p0_ut starts from 0 successes
  b_pa_x1    <-    b_pa[[n1]][1 + x1]                       # P(X1=x1)                                 , so   dbinom(x1    ,n1,pa)
  B_pa_s     <-  1-B_pa[[n2]][1 + (s-x1)]                   # P(X1+X2>=s+1)=1-P(X1+X2<=s)=1-P(X2<=s-x1), so 1-pbinom(s-x1-1,n2,pa)
  b_pa_s     <-    b_pa[[n2]][1 + (s-x1)]                   # P(X2=s)                                  , so   dbinom(s     ,n2,pa)
                                                            # B_pa_s + b_pa_s = P(X1+X2>=s). Note that calculating P(X1+X2>=s)=1-P(X1+X2<=s-1)
                                                            # is not possible here: if x1=s, then in P(X2<=s-x1-1), right side is negative, so no index
                                                            # in listings with binomial results
  alpha_temp <- B_p0_r1 + rev(cumsum(rev(b_p0_x1 * (b_p0_s+B_p0_s))))
  pi_temp    <- B_pa_r1 + rev(cumsum(rev(b_pa_x1 * (b_pa_s+B_pa_s))))
  list(N = n1 + n2, n1 = n1, n2 = n2, r1 = r_1, s = s, alpha_temp = alpha_temp, pi_temp = pi_temp)
}

#' @title The Sargent 2-stage function
#' @description This function calculates sample sizes of the Sargent 2-stage design.
#' @description The goal of a phase II trial is to make a preliminary determination regarding the activity and
#' tolerability of a new treatment and thus to determine whether the treatment warrants
#' further study in the phase III setting. \cr
#' This function calculates the sample size needed in a Sargent 2-stage design which is a
#' three-outcome design that allows for three outcomes: reject \eqn{H(0)}, reject \eqn{H(a)}, or reject neither.
#' @param p0 probability of the uninteresting response (null hypothesis \eqn{H0})
#' @param pa probability of the interesting response (alternative hypothesis Ha)
#' @param alpha Type I error rate \eqn{P(reject H0|H0)}
#' @param beta Type II error rate \eqn{P(reject Ha|Ha)}
#' @param eta \eqn{P(reject Ha|H0)}
#' @param pi \eqn{P(reject H0|Ha)}
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
#' \item r1: critical value for the first stage
#' \item r2: critical value for the second stage
#' \item eff: s/N
#' \item 90%CI_low: Result of call to getCI_Koyama. Confidence interval according to Koyama and Chen (1989)
#' \item 90%CI_high: Result of call to getCI_Koyama. Confidence interval according to Koyama and Chen (1989)
#' \item EN.p0: expected sample size under H0
#' \item PET.p0: probability of terminating the trial at the end of the first stage under H0
#' \item MIN: column indicating if the design is the minimal design
#' \item OPT: column indicating if the setting is the optimal design
#' \item ADMISS: column indicating if the setting is the admissible design
#' \item alpha: the actual alpha value which is smaller than \code{alpha_param + eps}
#' \item beta: the actual beta value where which is smaller than \code{beta_param + eps}
#' \item eta: the actual eta value which is smaller than \code{eta_param - eps}
#' \item pi: the actual pi value which is smaller than \code{pi_param - eps}
#' \item lambda:  1-(eta+alpha)
#' \item delta: 1-(beta+pi)
#' \item p0: your provided \code{p0} value
#' \item pa: your provided \code{pa} value
#' \item alpha_param: your provided \code{alpha} value
#' \item beta_param: your provided \code{beta} value
#' \item eta_param: your provided \code{eta} value
#' \item pi_param: your provided \code{pi} value
#' }
#' @details
#' if x1<=r1 --> stop futility \cr
#' if (x1+x2)<=r --> futility \cr
#' if (x1+x2)>=s --> efficacy \cr
#' @references Sargent DJ, Chan V, Goldberg RM. A three-outcome design for phase II clinical trials. Control Clin Trials. 2001;22(2):117-125. doi:10.1016/s0197-2456(00)00115-x
#' @export
#' @examples
#' samplesize <- sargent2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.1, eta = 0.8, pi = 0.8,
#'                             eps = 0.005, N_min = 15, N_max = 30)
#' plot(samplesize)
#'
#' \donttest{
#' data(data_sargent2)
#' test <- data_sargent2
#' samplesize <- sargent2stage(p0 = test$p0, pa = test$pa, alpha = test$alpha, beta = test$beta,
#'                             eta = test$eta, pi = test$pi,
#'                             eps = 0.005,
#'                             N_min = test$N_min, N_max = test$N_max)
#' optimal <- lapply(samplesize, FUN=function(x) subset(x, OPT == "Optimal"))
#' optimal <- data.table::rbindlist(optimal)
#' minimax <- lapply(samplesize, FUN=function(x) subset(x, MIN == "Minimax"))
#' minimax <- data.table::rbindlist(minimax)
#' }

sargent2stage <- function(p0, pa, alpha, beta, eta, pi, eps = 0, N_min, N_max, int=0, int_window=0.025){

  if(length(p0) > 1 && length(pa) > 1){
    results <- mapply(null = p0, alternative = pa, alpha = alpha, beta = beta, eta = eta, pi = pi, eps = eps, N_min = N_min, N_max = N_max, int = int, int_window = int_window,
                      FUN = function(null, alternative, alpha, beta, eta, pi, eps, N_min, N_max, int, int_window){
                        sargent2stage.default(p0 = null, pa = alternative, alpha = alpha, beta = beta, eta = eta, pi = pi, eps = eps, N_min = N_min, N_max = N_max, int = int, int_window = int_window)
                      }, SIMPLIFY = FALSE)
  }else{
    results <- sargent2stage.default(p0 = p0, pa = pa, alpha = alpha, beta = beta, eta = eta, pi = pi, eps = eps, N_min = N_min, N_max = N_max, int = int, int_window = int_window)
  }
}

sargent2stage.default <- function(p0, pa, alpha, beta, eta, pi, eps = 0, N_min, N_max, int=0, int_window=0.025) {

  # WARNING MESSAGES
  if (pa < p0) {
    stop("p0 should be smaller than pa")
  }

  if (N_min <1 ) {
    stop("N_min should be >=1")
  }

  EN.p0 <- EN.p0_N_min <- EN.p0_N_n1_min <- EN.p0_min <- EN.p0_N_min_int <- MIN <- N <- OPT <- NULL
  rowid <- n1 <- n2 <- r <- r1 <- r1min <- r2 <- s <- NULL

  #----------------------------------------------------------------------------------------------------------#
  # Get all possible scenarios for N, n1, r1, r2 and n2                                                      #
  # Note that this is the possible range of values: 0 <=r1<n1; r1+1<=r; r+2<=s<=n1+n2                        #
  #----------------------------------------------------------------------------------------------------------#

  # Get all N's for which there is a max r (1:N-2) for which P(X<=r|Ha)<=beta+eps, and select that max r
  # Note that this is not taking into account first stage. However, the probability of rejecting Ha
  # should be lower in the second stage, compared to a 1-stage design, as some cases already rejected at first stage,
  # meaning that an rmax, calculated using a 1-stage design, is the maximum
  #------------------------------------------------------------------------------------------------------------

  res0 <- lapply(N_min:N_max, FUN=function(a) cbind(N = a, rtemp = 1:(a - 2)))
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
  res1 <- mapply(a = res0$N, b = res0$rmax, FUN = function(a, b) cbind(N = a, n1 = (1:(a - 1)), rmax = b), SIMPLIFY = FALSE)
  res1 <- do.call(rbind, res1)
  res1 <- data.frame(res1)

  # Get for selected N's and n1's, r's (1:rmax)
  #--------------------------------------------
  res2 <- mapply(a = res1$N, b = res1$n1, c = res1$rmax,
                 FUN = function(a,b,c) cbind(N = a, n1 = b, rmax = c, r = (1:c)),
                 SIMPLIFY = FALSE)
  res2 <- do.call(rbind, res2)
  res2 <- data.frame(res2)
  res2$n2 <- res2$N - res2$n1

  # Calculate beta and eta
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
  res3 <- res2[, P_Sargent_reject_Ha(n1 = n1, n2 = n2, r = r, b_p0 = b_p0, B_p0 = B_p0, b_pa = b_pa, B_pa = B_pa), by = list(rowid)]
  res3 <- data.table::setDF(res3)
  res3$diff_beta <- res3$beta_temp - beta
  res4 <- res3[res3$diff_beta <= eps, ]
  res4$diff_eta <- eta - res4$eta_temp
  res5 <- res4[res4$diff_eta <= eps, ]
  res3 <- res5

  # Get for selected N's, n1's, r1's and r1's: s's ((r2+2):n)
  #----------------------------------------------------------
  res4 <- mapply(a = res3$N, b = res3$n1, c = res3$n2, d = res3$r1, e = res3$r2, f = res3$beta_temp, g = res3$eta_temp,
                 FUN = function(a, b, c, d, e, f, g) cbind(N = a, n1 = b, n2 = c, r1 = d, r2 = e, s = ((e + 2):a), beta_temp = f, eta_temp = g))
  res4 <- do.call(rbind, res4)
  res4 <- data.frame(res4)

  #----------------------------------------------------------#
  # Calculate pi and alpha for all scenarios (s needed)      #
  #----------------------------------------------------------#

  res4     <- data.table::setDT(res4)
  settings <- res4[, list(r1min = min(r1)), by = list(N, n1, n2, s)]
  settings$rowid <- seq_len(nrow(settings))
  res5 <- settings[, P_Sargent_reject_H0(n1 = n1, n2 = n2, r1min=r1min, s = s, b_p0 = b_p0, B_p0 = B_p0, b_pa = b_pa, B_pa = B_pa), by = list(rowid)]

  res5 <- merge(res4[, c("N", "n1", "n2", "r1", "r2", "s", "beta_temp", "eta_temp")], res5,
                  by = c("N", "n1", "n2", "r1", "s"), all.x = TRUE, all.y = FALSE, sort = FALSE)
  res5 <- data.table::setDF(res5)
  res5$diff_pi <- pi - res5$pi_temp
  res5 <- res5[res5$diff_pi <= eps, ]
  res5$diff_alpha <- res5$alpha_temp - alpha
  res5 <- res5[res5$diff_alpha <= eps, ]

  res5$r1<-as.numeric(res5$r1) # Make variable numeric

  if(nrow(res5) == 0){
    stop("No data satisfying the H0/Ha criteria")
  }
  res5$alpha <- alpha
  res5$beta  <- beta
  res5$eta   <- eta
  res5$pi    <- pi

  res5$PET.p0 <- pbinom(q = res5$r1, size = res5$n1, prob = p0, lower.tail = T)
  res5$EN.p0 <- res5$N - ((res5$N - res5$n1) * res5$PET.p0)

  res5 <- data.table::setDT(res5)

  res5 <- res5[, N_min := min(N)]
  res5 <- res5[, EN.p0_min := min(EN.p0)]
  res5 <- res5[, EN.p0_N_min := min(EN.p0), by = N]
  res5 <- res5[, EN.p0_N_n1_min := min(EN.p0) , by = list(N,n1)]

  if (int>0){
    res6_int <- res5[EN.p0 == EN.p0_N_n1_min & floor((int-int_window)*N)<=n1 & n1<=ceiling((int+int_window)*N), ]
    res6_int <- res6_int[, EN.p0_N_min_int := min(EN.p0) , by = list(N)]
    res6_int <- res6_int[EN.p0 == EN.p0_N_min_int   , ]
    res6_int$INTERIM=int;
  }
  res6 <- res5[EN.p0 == EN.p0_N_min, ]

  res6$OPT <- res6$MIN <- res6$ADMISS <- c("")
  res6$OPT[which(res6$EN.p0 == res6$EN.p0_min)] <- "Optimal"
  res6$MIN[which(res6$N == res6$N_min)] <- "Minimax" # Note: if multiple designs that meet the criteria for minimax:choose
  # design with minimal expected sample size under H0: "Optimal minimax design"

  res6 <- data.table::setDF(res6)

  # Get admissible designs
  y <- res6[, c("N", "EN.p0")]
  con.ind <- chull(y)[chull((y)) == cummin(chull((y)))]

  res <- res6[res6$N >= min(res6$N[res6$MIN == "Minimax"]) & res6$N <= max(res6$N[res6$OPT == "Optimal"]), ]
  res$ADMISS[which((rownames(res) %in% c(con.ind)) & (res$N > res$N_min) & (res$EN.p0 > res$EN.p0_min) & (res$N < unique(res$N[res$OPT == "Optimal"])))] <- "Admissible"

  if (int>0){
    res$EN_opt  <- "EN.p0_optimal"
    res6_int<-res6_int[res6_int$N >= min(res6$N[res6$MIN == "Minimax"]) & res6_int$N <= max(res6$N[res6$OPT == "Optimal"]), ]
    res<-merge(x = res, y = res6_int, by = c("N","n1","n2","r1","s","r2","beta_temp","eta_temp","rowid","alpha_temp","pi_temp","diff_pi","diff_alpha","alpha","beta","eta","pi",
                                             "PET.p0","EN.p0","N_min","EN.p0_min","EN.p0_N_min","EN.p0_N_n1_min"), all = TRUE)

    res$ADMISS[is.na(res$ADMISS)] <- "" # AFter merging: replace NA's with """
    res$MIN   [is.na(res$MIN)   ] <- ""
    res$OPT   [is.na(res$OPT)   ] <- ""
    res$EN_opt[is.na(res$EN_opt)] <- ""
  }

  # Calculate 1-2*alpha confidence interval, based on Koyama, Statistics in Medicine 2008

  res$eff <- paste0(res$s, "/", res$N, " (", 100 * round((res$s) / res$N, 3), "%)")
  CI <- mapply(a = res$s, b = res$r1, c = res$n1, d = res$N,
               FUN = function(a, b, c, d) getCI_Koyama(k = a, r1 = b, n1 = c, n = d, alpha = alpha, precision = 4))
  res$CI_low  <- 100 * unlist(CI[rownames(CI) == "CI_low", ])
  res$CI_high <- 100 * unlist(CI[rownames(CI) == "CI_high", ])
  res <- data.table::setnames(res,
                              old = c("alpha", "beta", "eta", "pi"),
                              new = c("alpha_param", "beta_param", "eta_param", "pi_param"))
  res <- data.table::setnames(res,
                              old = c("alpha_temp", "beta_temp", "eta_temp", "pi_temp"),
                              new = c("alpha", "beta", "eta", "pi"))
  res$lambda <- 1 - (res$eta + res$alpha)
  res$delta  <- 1 - (res$beta + res$pi)
  res <- data.table::setDF(res)

  if (int>0){
	  res <- cbind(design_nr=1:dim(res)[1],
	               res[, c("r1", "n1", "r2", "s", "n2", "N", "eff", "CI_low", "CI_high", "EN.p0", "PET.p0", "MIN", "OPT", "ADMISS", "EN_opt","INTERIM", "alpha", "beta", "eta", "pi", "lambda", "delta")],
	               p0 = p0, pa = pa,
	               res[, c("alpha_param", "beta_param", "eta_param", "pi_param")])
	  }

  if (int==0){
	  res <- cbind(design_nr=1:dim(res)[1],
	               res[, c("r1", "n1", "r2", "s", "n2", "N", "eff", "CI_low", "CI_high", "EN.p0", "PET.p0", "MIN", "OPT", "ADMISS", "alpha", "beta", "eta", "pi", "lambda", "delta")],
	               p0 = p0, pa = pa,
	               res[, c("alpha_param", "beta_param", "eta_param", "pi_param")])
	  }

  res <- data.table::setnames(res,
                              old = c("CI_low", "CI_high"),
                              new = c(paste0(100 - 2 * 100 * alpha, "%CI_low"), paste0(100 - 2 * 100 * alpha, "%CI_high")))
  attr(res, "inputs") <- list(p0 = p0, pa = pa, alpha = alpha, beta = beta, eta = eta, pi = pi, eps = eps, N_min = N_min, N_max = N_max)
  class(res) <- c("2stage", "sargent", "data.frame")
  res
}

# TEST (Sargent DJ, Goldberg RM. A Three-Outcome Design for Phase II Clinical Trials. Controlled Clinical Trials 22:117-125
#--------------------------------------------------------------------------------------------------------------------------
# test_0_a<- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,2),pa=rep((a+0.15),2),alpha=c(0.1,0.05),beta=c(0.1,0.1),eta=c(0.8,0.8),pi=c(0.8,0.8)),a=c(0.05,0.1,0.2,0.3,0.4,0.5),SIMPLIFY=F)))
# test_0_b<- data.frame(do.call("rbind", mapply(function(a,b) cbind(p0=rep(a,2),pa=rep((a+0.2 ),2),alpha=c(0.1,0.05),beta=c(0.1,0.1),eta=c(0.8,0.8),pi=c(0.8,0.8)),a=c(0.05,0.1,0.2,0.3,0.4,0.5),SIMPLIFY=F)))
# test_0<-cbind(rbind(test_0_a,test_0_b),data.frame(N_min=c(23,23,26,31,37,47,44,58,48,64,46,62,12,12,17,17,20,28,24,32,24,33,24,33),N_max=c(31,31,37,45,52,58,58,69,56,74,56,72,21,21,28,30,29,38,32,41,36,41,37,43)))
#
# for (i in 1:dim(test_0)[1]){
#   res<- PhIIdesign::sargent2stage(p0=test_0[i,]$p0,pa=test_0[i,]$pa,alpha=test_0[i,]$alpha,beta=test_0[i,]$beta,eta=test_0[i,]$eta,pi=test_0[i,]$pi,eps=0.005,N_min=test_0[i,]$N_min,N_max=test_0[i,]$N_max)
#
#   if (i==1) {test_list_O      <-list(res[OPT=="Optimal",])}   # Create list with first dataset
#   if (i!=1) {test_list_O[[i]] <-res[OPT=="Optimal",] }        # Next iterations: append dataset
#
#   if (i==1) {test_list_M      <-list(res[MIN=="Minimax",])}   # Create list with first dataset
#   if (i!=1) {test_list_M[[i]] <-res[MIN=="Minimax",] }        # Next iterations: append dataset
#
# }
# test_O<- data.frame(do.call("rbind",test_list_O))
# test_O<-test_O[,c("p0","pa","alpha_param","beta_param","eta_param","pi_param","alpha","beta","eta","pi","n1","n2","r1","r2","s","N","EN.p0")]
#
# test_M<- data.frame(do.call("rbind",test_list_M))
# test_M<-test_M[,c("p0","pa","alpha_param","beta_param","eta_param","pi_param","alpha","beta","eta","pi","n1","n2","r1","r2","s","N","EN.p0")]
