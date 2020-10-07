
#' @title Plot function for 1- and 2-stage Sargent design
#' @description Plots the probability of rejecting \eqn{H0} ('Positive result') and rejecting Ha \cr
#' ('Negative result') for different true response binomial parameters \cr
#' if the input is an object returned by \code{\link{sargent2stage}}, interactively choose \cr
#' the design number
#' @param x an object returned by \code{\link{sargent1stage}} or \code{\link{sargent2stage}}
#' @param design_nr design_nr returned by \code{\link{sargent1stage}} or \code{\link{sargent2stage}} \cr
#' default is 1 (for \code{\link{sargent1stage}}, only one design is outputted)
#' @param main title of the graph, passed on to \code{plot}
#' @param xlab x-axis label of the graph, passed on to \code{plot}
#' @param ylab y-axis label of the graph, passed on to \code{plot}
#' @param ... other arguments passed on to \code{plot}
#' @export
#' @examples
#' \donttest{
#' result_1<- sargent1stage(p0 = 0.2, pa = 0.35, alpha = 0.1, beta = 0.1, eta = 0.8, pi = 0.8,
#'                          eps = 0.005, N_min = 35, N_max = 50)
#' plotsargent(result_1)
#'
#' result_2 <- sargent2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.1,
#'                           eta = 0.8, pi = 0.8,
#'                           eps = 0.005, N_min = 15, N_max = 30)
#' plotsargent(result_2)
#' }

plotsargent <- function(x,
                        design_nr=1,
                        main = "Sargent's 3-outcome 1-stage design",
                        xlab = "True response probability",
                        ylab, ...) {

  x<-x[x$design_nr==design_nr,]

  #res: truep, rejectHA,nodecision,reject H0
  Res=NULL
  ptrue=seq(max(0,x$p0-0.1),min(x$pa+0.2,1),by=0.01)
  for (i in 1:length(ptrue)){
    p=ptrue[i]

    if ("r1" %in% names(x)){
      rejectH0=prob_reject_H0(n1=x$n1,n2=x$n2,r1=x$r1,s =x$s ,p=p)
      rejectHA=prob_reject_Ha(n1=x$n1,n2=x$n2,r1=x$r1,r2=x$r2,p=p)
    }

    if (!"r1" %in% names(x)){
      rejectH0=pbinom(x$s-1, x$N, p, lower.tail = FALSE)
      rejectHA=pbinom(x$r,x$N, p)
    }

    Res_p1=c(p, rejectHA, 1-(rejectH0+rejectHA),rejectH0)
    Res=rbind(Res,Res_p1)
  }

  colnames(Res)<-c("Truep","ProbrejHA","Probnodec","ProbrejH0")
  rownames(Res)=NULL
  plotdata=as.data.frame(Res)
  cols=c("green2","green4","red2","red4")
  par(mar = c(5, 4, 4, 4) + 0.3)
  #should have symbols instead of <= and >=
  plot(plotdata$Truep,plotdata$ProbrejHA,type="l",xlim=c(x$p0,min(x$pa+0.2)),ylim=c(0,1),
       xlab=xlab,ylab="Probability negative result",
       main=paste0(main,"\n Parameters: H0:p<=",x$p0,", Ha:p>=",x$pa,", alpha=",x$alpha_param,", beta=",x$beta_param,", eta=",x$eta_param, ", pi=",x$pi_param,"\n Design N=",x$N,", r=",x$r," and s=",x$s),yaxt="n" ,col=cols[4])
  axis(side = 2, at = pretty(range(plotdata$ProbrejHA+plotdata$Probnodec)),las=1,
       labels=sort(pretty(range(plotdata$ProbrejHA+plotdata$Probnodec)),decreasing=F))
  par(new = TRUE)
  plot(plotdata$Truep,plotdata$ProbrejHA+plotdata$Probnodec, type="l",
       axes = FALSE, xlab = "", ylab = "",xlim=c(x$p0,min(x$pa+0.2)),ylim=c(0,1),col=cols[2])
  axis(side = 4, at = pretty(range(plotdata$ProbrejHA+plotdata$Probnodec)),las=1,
       labels=sort(pretty(range(plotdata$ProbrejHA+plotdata$Probnodec)),decreasing=T))
  p <- par('usr')
  text(p[2]+0.04, mean(p[3:4]), labels = 'Probability positive results', xpd = T, srt = -90)

  polygon(c(plotdata$Truep[1],plotdata$Truep),c(0,plotdata$ProbrejHA),
          border = NA,col=cols[4])
  polygon(c(plotdata$Truep,rev(plotdata$Truep)),
          c((plotdata$ProbrejHA+plotdata$Probnodec),rep(1,length(plotdata$Truep))),
          border=NA, col=cols[2])
  # grid()
  respo=round(subset(plotdata,as.numeric(plotdata$Truep)==x$p0)[2:4]*100,1)
  respa=round(subset(plotdata,as.numeric(plotdata$Truep)==x$pa)[2:4]*100,1)
  #I don't understand why the above doesn't work, it works for p0 but not for pa
  #Can you solve this Jan?
  respa=round(plotdata[36,2:4]*100,1)
  text(x=x$p0 ,y=0.08,pos=4, "Negative result",font=2,col=cols[3],cex=1.5)
  text(x=x$p0,y=0.33,pos=4,paste0("Under H0 (p=",x$p0,"):"),col=cols[3],cex=1.1)
  text(x=x$p0,y=0.30,pos=4,paste0(respo[1],"% Prob negative"),col=cols[3],cex=1.1)
  text(x=x$p0,y=0.27,pos=4,paste0(respo[2],"% Prob inconclusive"),col=cols[3],cex=1.1)
  text(x=x$p0,y=0.24,pos=4,paste0(respo[3],"% Prob positive (alpha)"),col=cols[3],cex=1.1)
  text(x=x$pa+0.08,y=0.95, pos=4,"Positive result",font=2,col=cols[1],cex=1.5)
  text(x=x$pa,y=0.70,pos=4,paste0("Under Ha (p=",x$pa,"):"),col=cols[1],cex=1.1)
  text(x=x$pa,y=0.67,pos=4,paste0(respa[1],"% Prob negative"),col=cols[1],cex=1.1)
  text(x=x$pa,y=0.64,pos=4,paste0(respa[2],"% Prob inconclusive"),col=cols[1],cex=1.1)
  text(x=x$pa,y=0.61,pos=4,paste0(respa[3],"% Prob positive (power)"),col=cols[1],cex=1.1)
  abline(v=c(x$p0,x$pa),lty=2,col=cols[c(3,1)])

  return(as.data.frame(Res))
}


