

P_Fleming2st_reject_Ha <- function(n1, n2, r, a, b_p0, B_p0, b_pa, B_pa){
  ## r is the total number of seen cases on both n1 and n2 together
  ## Calculate alpha/beta for all values up to r (see formula Section 2.2, paper Mander and Thompson)

  r_1 <- (max(0,r-n2)):(max(0,min(a-2, r-1))) # for min: - there should be x's for which r1<x1<a, so minimum for a should be 2
                                              #          - if r>n2, and knowing that you must have (r+1) successes to reject H0 after stage 2,
                                              #            then you must already have >= (r+1)-n2 successes in stage 1 with continuation, so r1 should be minimal
                                              #            (r-n2), and with one x1 between r1 and a, the minimal value for a is thus (r-n2)+2
                                              # for max: r1<=a-2, as there must be at least one value x1, for which r1<x1<a, so trial can continue to stage 2
                                              #          also logical that r1<=r-1, as you need r1+1 to proceed to second stage.
                                              #          if r1=r than in second stage X1+X2>r, so second stage not needed anymore

  x1  <- r_1 + 1  # Function is to calculate P(reject Ha), so going from (r_1 + 1) to min(n1,r)
  #   'n1' is logical. 'r': if X1>r then you cannot reject Ha. [see also previous step]

  # calculate alpha for each (r1,x1)
  B_p0_r1    <- B_p0[[n1]][1 + r_1]     # P(X1<=r1|n1,p0)   so pbinom(r_1 ,n1,p0), indexing with "+1", as B_p0 starts from 0 successes
  b_p0_x1    <- b_p0[[n1]][1 + x1]      # P(X1=x1|n1,p0)    so dbinom(x1  ,n1,p0), indexing with "+1", as b_p0 starts from 0 successes
  B_p0_r2    <- B_p0[[n2]][1 + (r-x1)]  # P(X2<=r-x1|n2,p0) so pbinom(r-x1,n2,p0), indexing with "+1", as B_p0 starts from 0 successes
  eta_temp   <- B_p0_r1 + rev(cumsum(rev(b_p0_x1 * B_p0_r2))) # e.g. for r=2: r_1=0: P(X1<=0) + P(X1=1)P(X2<=2-1) + P(X1=2)P(X2<=2-2)
  #               r_1=1: P(X1<=1) +                     P(X1=2)P(X2<=2-2)
  # calculate beta for each (r1,x1)
  B_pa_r1    <- B_pa[[n1]][1 + r_1]
  b_pa_x1    <- b_pa[[n1]][1 + x1]
  B_pa_r2    <- B_pa[[n2]][1 + (r-x1)]
  beta_temp  <- B_pa_r1 + rev(cumsum(rev(b_pa_x1 * B_pa_r2)))
  list(N = n1 + n2, n1 = n1, r1 = r_1, a=a, n2 = n2, r2 = r, alpha_temp = 1 - eta_temp, beta_temp = beta_temp)
}

#' @title Fleming 2-stage function
#' @description
#' The primary objective of a phase II clinical trial of a new drug or regimen is to determine whether it has sufficient biological activity
#' against the disease under study to warrant more extensive development. \cr
#' This function calculates the sample size needed in a Fleming 2-stage design which is a
#' two-stage design that is optimal in the sense that the expected sample size is minimized if the
#' regimen has low activity subject to constraints upon the size of the type 1 and type 2 errors. \cr
#' Two-stage designs which minimize the maximum sample size are also determined. \cr
#' This type of design also allows stopping for efficacy
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
#' \item a: ("a" for "acceptance") threshold for "accepting" Ha= stop for efficacy at first stage
#' \item eff: (r2 + 1)/N
#' \item 90%CI_low: exact 1-2*alpha confidence interval lower limit
#' \item 90%CI_high: exact 1-2*alpha confidence interval upper limit
#' \item EN.p0: expected sample size under H0
#' \item PET.p0: probability of terminating the trial for futility at the end of the first stage under H0
#' \item EN.pa: expected sample size under Ha
#' \item PET.pa: probability of terminating the trial for efficacy at the end of the first stage under Ha
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
#' if x1<=r1 --> stop futility at first stage \cr
#' if r1<x1<a --> proceed to second dtage \cr
#' if x1>=a --> stop efficacy at first stage \cr
#' if (x1+x2)<=r2 --> futility at second stage \cr
#' if (x1+x2)> r2 --> efficacy at second stage \cr
#' So PET.H0=P_H0(X1<=r1)+P_H0(X1>=a1)
#' @references Mander AP, Thompson SG. Two-stage designs optimal under the alternative hypothesis for phase II cancer clinical trials.
#'     Contemporary Clinical Trials 2010;31:572–578
#'     Qin F et al. Optimal, minimax and admissible two-stage design for phase II oncology clinical trials. BMC Medical Research Methodology 2020;20:126
#' @export
#' @examples
#' samplesize <- fleming2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.2,
#'                           eps = 0.005, N_min = 1, N_max = 50)

fleming2stage <- function(p0, pa, alpha, beta, eps = 0, N_min, N_max, int=0, int_window=0.025){

  if(length(p0) > 1 && length(pa) > 1){
    results <- mapply(null = p0, alternative = pa, alpha = alpha, beta = beta, eps = eps, N_min = N_min, N_max = N_max, int = int, int_window = int_window,
                      FUN = function(null, alternative, alpha, beta, eps, N_min, N_max, int, int_window){
                        fleming2stage.default(p0 = null, pa = alternative, alpha = alpha, beta = beta, eps = eps, N_min = N_min, N_max = N_max, int = int, int_window = int_window)
                      }, SIMPLIFY = FALSE)
  }else{
    results <- fleming2stage.default(p0 = p0, pa = pa, alpha = alpha, beta = beta, eps = eps, N_min = N_min, N_max = N_max, int = int, int_window = int_window)
  }
}

# p0 = 0.1; pa = 0.3; alpha = 0.05; beta = 0.2; eps = 0.005; N_min = 20; N_max = 21

fleming2stage.default <- function(p0, pa, alpha, beta, eps = 0, N_min, N_max, int=0, int_window=0.025) {

  # WARNING MESSAGES
  if (pa < p0) {
    stop("p0 should be smaller than pa")
  }

  if (N_min <1 ) {
    stop("N_min should be >=1")
  }

  EN.p0 <- EN.p0_N_min <- EN.p0_N_n1_min <- EN.p0_min <- EN.p0_N_min_int <- EN.pa <- MIN <- N <- OPT <- NULL
  rowid <- n1 <- n2 <- r <- r2 <- a <- NULL

  #----------------------------------------------------------------------------------------------------------#
  # Get all possible scenarios for N, n1, r1, r2 and n2                                                      #
  # Note that this is the possible range of values: 0 <=r1<n1; r1+1<=r; r<=n1+n2                             #
  #----------------------------------------------------------------------------------------------------------#

  # Get all N's for which there is a max r (1:N) for which P(X<=r|Ha)<=beta+eps, and select that max r
  # Note that this is not taking into account first stage. However, the probability of rejecting Ha
  # should be lower in the second stage, compared to a 1-stage design, as some cases already rejected at first stage,
  # meaning that an rmax, calculated using a 1-stage design, is a good maximum (the higher r, the higher the probability to reject Ha)
  #-----------------------------------------------------------------------------------------------------------------------------------

  res0 <- lapply(N_min:N_max, FUN=function(N) cbind(N = N, rtemp = 1:N))
  res0 <- do.call(rbind, res0)
  res0 <- data.frame(res0)
  res0$betamax <- pbinom(q = res0$rtemp, size = res0$N, prob = pa, lower.tail = TRUE)
  res0 <- res0[!is.na(res0$betamax) & res0$betamax <= (beta + eps), ]
  # Select all possible N"s: there needs to be
  # at least one r with P(X<=r|Ha)<=beta+eps
  res0 <- aggregate(res0$rtemp, by = list(res0$N), max)
  names(res0) <- c("N", "rmax")

  # Get for selected N's all possible n1's (2:N-1)  # there should be x's for which r1<x1<a, so minimum for a should be 2, and thus minimum for n1 should be 2
  #------------------------------------------------------------------------------------------------------------------------------------------------------------
  res1 <- mapply(N = res0$N, rmax = res0$rmax, FUN = function(N, rmax) cbind(N = N, n1 = (2:(N - 1)), rmax = rmax), SIMPLIFY = FALSE)
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


  # Get for selected N's, n1's, r's: a's (max(2,r-n2+2)):n1)  Note: for min: - there should be x's for which r1<x1<a, so minimum for a should be 2
  #                                                                          - if r>n2, and knowing that you must have (r+1) successes to reject H0 after stage 2,
  #                                                                            then you must already have >= (r+1)-n2 successes in stage 1 with continuation, so r1 should be minimal
  #                                                                            (r-n2), and with one x1 between r1 and a, the minimal value for a is thus (r-n2)+2
  #                                                                for max: a<=n1 that's logical
  #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  res3 <- mapply(N = res2$N, n1 = res2$n1, rmax = res2$rmax,r=res2$r, n2=res2$n2,
                 FUN = function(N, n1, rmax, r, n2) cbind(N = N, n1 = n1, rmax = rmax, r = r, n2 = n2, a=(max(2,r-n2+2)):n1),
                 SIMPLIFY = FALSE)
  res3 <- do.call(rbind, res3)
  res2 <- data.frame(res3)


  # Calculate beta and alpha
  #-------------------------

  ###   use a lookup table for formula Section 2.2 from paper Mander and Thompson
  ###   b: probability mass
  ###   B: cumulative probability
  nmax <- N_max
  b_p0 <- lapply(1:nmax, FUN = function(n) dbinom(0:nmax, size = n, prob = p0))
  b_pa <- lapply(1:nmax, FUN = function(n) dbinom(0:nmax, size = n, prob = pa))
  B_p0 <- lapply(1:nmax, FUN = function(n) pbinom(0:nmax, size = n, prob = p0, lower.tail = TRUE))
  B_pa <- lapply(1:nmax, FUN = function(n) pbinom(0:nmax, size = n, prob = pa, lower.tail = TRUE))

  res2 <- data.table::setDT(res2)
  res2$rowid <- seq_len(nrow(res2))
  res3 <- res2[, P_Fleming2st_reject_Ha(n1 = n1, n2 = n2, r = r, a=a, b_p0 = b_p0, B_p0 = B_p0, b_pa = b_pa, B_pa = B_pa), by = list(rowid)]
  res3 <- data.table::setDF(res3)
  res3$diff_beta <- res3$beta_temp - beta
  res4 <- res3[res3$diff_beta <= eps, ]
  res4$diff_alpha <- res4$alpha_temp - alpha
  res5 <- res4[res4$diff_alpha <= eps, ]

  res5[2:7] <- lapply(res5[2:7], as.numeric) # convert integer into numeric

  if(nrow(res5) == 0){
    stop("No data satisfying the H0/Ha criteria")
  }

  res5$alpha <- alpha
  res5$beta <- beta

  res5$PET.p0 <- pbinom(q = res5$r1, size = res5$n1, prob = p0, lower.tail = TRUE) + pbinom(q = res5$a-1, size = res5$n1, prob = p0, lower.tail = FALSE)
  res5$PET.pa <- pbinom(q = res5$r1, size = res5$n1, prob = pa, lower.tail = TRUE) + pbinom(q = res5$a-1, size = res5$n1, prob = pa, lower.tail = FALSE)

  res5$EN.p0  <- res5$N - ((res5$N - res5$n1) * res5$PET.p0)
  res5$EN.pa  <- res5$N - ((res5$N - res5$n1) * res5$PET.pa)

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
    res<-merge(x = res, y = res6_int, by = c("N","n1","n2","r1","a","r2","alpha_temp","beta_temp","diff_beta","diff_alpha","alpha","beta","PET.p0","EN.p0","PET.pa","EN.pa","N_min","EN.p0_min","EN.p0_N_min","EN.p0_N_n1_min","rowid"), all = TRUE)

    res$ADMISS[is.na(res$ADMISS)] <- "" # AFter merging: replace NA's with """
    res$MIN   [is.na(res$MIN)   ] <- ""
    res$OPT   [is.na(res$OPT)   ] <- ""
    res$EN_opt[is.na(res$EN_opt)] <- ""
  }


  # Calculate exact 1-2*alpha confidence interval

  res$eff <- paste0(res$r2 + 1, "/", res$N, " (", 100 * round((res$r2 + 1) / res$N, 3), "%)")
  CI <- mapply(a = res$r2 + 1, b = res$N, FUN = function(a, b) binom::binom.confint(a,b,conf.level=1-2*alpha,methods="exact"))
  res$CI_low  <- round(100 * unlist(CI[rownames(CI) == "lower", ]),2)
  res$CI_high <- round(100 * unlist(CI[rownames(CI) == "upper", ]),2)

  res <- data.table::setnames(res,
                              old = c("alpha", "beta"),
                              new = c("alpha_param", "beta_param"))
  res <- data.table::setnames(res,
                              old = c("alpha_temp", "beta_temp"),
                              new = c("alpha", "beta"))

  if (int>0){
  res <- cbind(design_nr=1:dim(res)[1],
               res[, c("r1", "a", "n1", "r2", "n2", "N", "eff", "CI_low", "CI_high", "EN.p0", "PET.p0", "EN.pa", "PET.pa", "MIN", "OPT", "ADMISS","EN_opt","INTERIM", "alpha", "beta")],
               p0 = p0, pa = pa,
               res[, c("alpha_param", "beta_param")])
  }
  if (int==0){
    res <- cbind(design_nr=1:dim(res)[1],
                 res[, c("r1", "a", "n1", "r2", "n2", "N", "eff", "CI_low", "CI_high", "EN.p0", "PET.p0", "EN.pa", "PET.pa", "MIN", "OPT", "ADMISS","alpha", "beta")],
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
