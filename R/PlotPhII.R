
#' @title Plot function for Probability of conclusions of a PhII design
#' @description Plots the probability of rejecting \eqn{H0} ('Positive result') and rejecting Ha \cr
#' ('Negative result') for different true response binomial parameters \cr
#' if the input is an object returned by \code{\link{fleming1stage}}, \code{\link{simon2stage}},
#' \code{\link{simon2stage}}, \code{\link{sargent1stage}} or \code{\link{sargent2stage}}, interactively choose \cr
#' the design number
#' @param x an object returned by \code{\link{fleming1stage}}, \code{\link{simon2stage}},
#' \code{\link{simon2stage}}, \code{\link{sargent1stage}} or \code{\link{sargent2stage}} \cr
#' @param design_nr design_nr returned by \code{\link{fleming1stage}}, \code{\link{simon2stage}},
#' \code{\link{simon2stage}}, \code{\link{sargent1stage}} or \code{\link{sargent2stage}} \cr
#' default is 1 (for \code{\link{sargent1stage}}, only one design is outputted)
#' @param main title of the graph, passed on to \code{plot}
#' @param xlab x-axis label of the graph, passed on to \code{plot}
#' @param ylab y-axis label of the graph, passed on to \code{plot}
#' @param grid add grid to plot ("Y" or "N")
#' @param ... other arguments passed on to \code{plot}
#' @references Sargent DJ, Chan V, Goldberg RM. A three-outcome design for phase II clinical \cr
#' trials. Control Clin Trials. 2001;22(2):117-125. doi:10.1016/s0197-2456(00)00115-x \cr
#' Simon R. Optimal two-Stage Designs for Phase II Clinical Trials. \cr
#' Control Clin Trials. 1989;10:1-10
#' @export
#' @examples
#' \donttest{
#' fleming1 <- fleming1stage(p0 = 0.45, pa = 0.7, alpha = 0.05, beta = 0.2)
#' plotPhII(fleming1)
#' simon2<- simon2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.2,
#'                      eps = 0.005, N_min = 1, N_max = 50)
#' plotPhII(simon2)
#'
#' sargent1 <- sargent1stage(p0 = 0.2, pa = 0.35, alpha = 0.1, beta = 0.1, eta = 0.8, pi = 0.8,
#'                          eps = 0.005, N_min = 35, N_max = 50)
#' plotPhII(sargent1)
#'
#' sargent2 <- sargent2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.1,
#'                           eta = 0.8, pi = 0.8,
#'                           eps = 0.005, N_min = 15, N_max = 30)
#' plotPhII(sargent2,design_nr=1)
#' }

plotPhII <- function(x,design_nr=1, main = "PhII design", xlab = "True response probability", grid="Y",ylab, ...) {

  x<-x[x$design_nr==design_nr,]

  ptrue<-seq(max(0,x$p0-0.1),min(x$pa+0.2,1),by=0.01)
  ProbrejHA<-Probnodec<-ProbrejH0<-rep(NA,length(ptrue))

  for (i in 1:length(ptrue)){

    p=ptrue[i]

    if ("r1" %in% names(x) & "s" %in% names(x)){  # for Sargent 2-stage
      ProbrejH0[i]=sargent2stage_prob_reject_H0(n1=x$n1,n2=x$n2,r1=x$r1,s =x$s ,p=p)
      ProbrejHA[i]=sargent2stage_prob_reject_Ha(n1=x$n1,n2=x$n2,r1=x$r1,r2=x$r2,p=p)
    }

    if (!"r1" %in% names(x) & "s" %in% names(x)){  # for Sargent 1-stage
      ProbrejH0[i]=pbinom(x$s-1, x$N, p, lower.tail = FALSE)
      ProbrejHA[i]=pbinom(x$r,x$N, p)
    }

    if (!"r1" %in% names(x) & !"s" %in% names(x)){  # for Fleming 1-stage
      ProbrejHA[i]=pbinom(x$r,x$N, p)
      ProbrejH0[i]=1-ProbrejHA[i]
    }

    if ("r1" %in% names(x) & !"s" %in% names(x)){  # For Simon 2-stage
      ProbrejHA[i]=simon2stage_prob_reject_Ha(n1=x$n1,n2=x$n2,r1=x$r1,r=x$r2,p=p)
      ProbrejH0[i]=1-ProbrejHA[i]
    }

    Probnodec[i]<-1-(ProbrejH0[i]+ProbrejHA[i])

  }

  cols=c("green2","green4","red2","red4")
  par(mar = c(5, 4, 4, 4) + 0.3)

  if ("r1" %in% names(x) & "s" %in% names(x)){  # for Sargent 2-stage
    plot(ptrue,ProbrejHA,type="l",ylim=c(0,1),xlab=xlab,ylab="Probability negative result",
         main=paste0(main,"\n Parameters: H0:p\u2264",x$p0,", Ha:p\u2265",x$pa,
                     ", alpha=",x$alpha_param,
                     ", beta=" ,x$beta_param,
                     ", eta="  ,x$eta_param,
                     ", pi="   ,x$pi_param,
                     "\n Design N=",x$N,", n1=",x$n1,", r1=",x$r1,", r2=",x$r2," and s=",x$s),yaxt="n" ,col=cols[4])
  }

  if (!"r1" %in% names(x) & "s" %in% names(x)){  # for Sargent 1-stage
    plot(ptrue,ProbrejHA,type="l",ylim=c(0,1),xlab=xlab,ylab="Probability negative result",
         main=paste0(main,"\n Parameters: H0:p\u2264",x$p0,", Ha:p\u2265",x$pa,
                     ", alpha=",x$alpha_param,
                     ", beta=" ,x$beta_param,
                     ", eta="  ,x$eta_param,
                     ", pi="   ,x$pi_param,
                     "\n Design N=",x$N,", r=",x$r," and s=",x$s),yaxt="n" ,col=cols[4])
  }

  if (!"r1" %in% names(x) & !"s" %in% names(x)){  # for Fleming 1-stage
    plot(ptrue,ProbrejHA,type="l",ylim=c(0,1),xlab=xlab,ylab="Probability negative result",
         main=paste0(main,"\n Parameters: H0:p\u2264",x$p0,", Ha:p\u2265",x$pa,
                     ", alpha=",x$alpha_param,
                     ", beta=" ,x$beta_param,
                     "\n Design N=",x$N,", r=",x$r),yaxt="n" ,col=cols[4])
  }

  if ("r1" %in% names(x) & !"s" %in% names(x)){  # For Simon 2-stage
    plot(ptrue,ProbrejHA,type="l",ylim=c(0,1),xlab=xlab,ylab="Probability negative result",
         main=paste0(main,"\n Parameters: H0:p\u2264",x$p0,", Ha:p\u2265",x$pa,
                     ", alpha=",x$alpha_param,
                     ", beta=" ,x$beta_param,
                     "\n Design N=",x$N,", n1=",x$n1,", r1=",x$r1,", r2=",x$r2),yaxt="n" ,col=cols[4])
  }

  axis(side = 2, at = seq(0,1,0.1),las=1,labels=seq(0,1,0.1))
  par(new = TRUE)
  plot(ptrue,ProbrejHA+Probnodec, type="l", axes = FALSE, xlab = "", ylab = "",ylim=c(0,1),col=cols[2])
  axis(side = 4, at = seq(0,1,0.1),las=1,labels=seq(1,0,-0.1))
  mtext("Probability positive results", side = 4, line = 3)

  polygon(c(ptrue[1],ptrue,ptrue[length(ptrue)]),c(0,ProbrejHA,ProbrejHA[length(ptrue)]), border = NA,col=cols[4])
  polygon(c(ptrue,rev(ptrue)),c(ProbrejHA+Probnodec,rep(1,length(ptrue))), border=NA, col=cols[2])

  if (grid=="Y"){grid()}

  text(x=x$p0 ,y=0.08,pos=4, "Negative result",font=2,col=cols[3],cex=1.5)
  text(x=x$p0,y=0.33,pos=4,paste0("Under H0 (p=",x$p0,"):"),col=cols[3],cex=1.1)
  text(x=x$p0,y=0.30,pos=4,paste0(round(100*ProbrejHA[as.character(ptrue)==as.character(x$p0)],1),"% Prob negative (eta)"),col=cols[3],cex=1.1)
  text(x=x$p0,y=0.27,pos=4,paste0(round(100*Probnodec[as.character(ptrue)==as.character(x$p0)],1),"% Prob inconclusive"),col=cols[3],cex=1.1)
  text(x=x$p0,y=0.24,pos=4,paste0(round(100*ProbrejH0[as.character(ptrue)==as.character(x$p0)],1),"% Prob positive (alpha)"),col=cols[3],cex=1.1)
  text(x=x$pa+0.08,y=0.95, pos=4,"Positive result",font=2,col=cols[1],cex=1.5)
  text(x=x$pa,y=0.70,pos=4,paste0("Under Ha (p=",x$pa,"):"),col=cols[1],cex=1.1)
  text(x=x$pa,y=0.67,pos=4,paste0(round(100*ProbrejHA[as.character(ptrue)==as.character(x$pa)],1),"% Prob negative (beta)"),col=cols[1],cex=1.1)
  text(x=x$pa,y=0.64,pos=4,paste0(round(100*Probnodec[as.character(ptrue)==as.character(x$pa)],1),"% Prob inconclusive"),col=cols[1],cex=1.1)
  text(x=x$pa,y=0.61,pos=4,paste0(round(100*ProbrejH0[as.character(ptrue)==as.character(x$pa)],1),"% Prob positive (power)"),col=cols[1],cex=1.1)
  abline(v=c(x$p0,x$pa),lty=2,col=cols[c(3,1)])

  return(data.frame(ptrue=ptrue,ProbrejHA=ProbrejHA,Probnodec=Probnodec,ProbrejH0=ProbrejH0))
}


