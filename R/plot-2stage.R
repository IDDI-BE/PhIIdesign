#' @title Plot function for 2 stage sample size calculation
#' @description Plots the expected sample size under H0 versus the maximum sample size
#' @param x an object returned by \code{\link{simon2stage}} or \code{\link{sargent2stage}}
#' @param main title of the graph, passed on to \code{plot}
#' @param xlab x-axis label of the graph, passed on to \code{plot}
#' @param ylab y-axis label of the graph, passed on to \code{plot}
#' @param ... other arguments passed on to \code{plot}
#' @export
#' @examples
#' samplesize <- simon2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.2,
#'                           eps = 0.005, N_min = 1, N_max = 50)
#' plot(samplesize)
#'
#' samplesize <- sargent2stage(p0 = 0.1, pa = 0.3, alpha = 0.05, beta = 0.1,
#'                             eta = 0.8, pi = 0.8,
#'                             eps = 0.005, N_min = 15, N_max = 30)
#' plot(samplesize)
plot.2stage <- function(x,
                        main = "Two-stage Designs",
                        xlab = "Maximum Sample Size N",
                        ylab, ...){
  if(missing(ylab)){
    ylab <- paste("E( N | ", attr(x, "inputs")$p0, " )", sep = "")
    ylab <- force(ylab)
    ylab <- as.expression(ylab)
  }

  if(is.null(x$EN_opt)){
    plot(x$N, x$EN.p0, type = "l", lwd=2, xlab = xlab, ylab = ylab, main = main, ...)
  }

  if(!is.null(x$EN_opt)){
    plot(x[which(x$EN_opt=="EN.p0_optimal"),]$N, x[which(x$EN_opt=="EN.p0_optimal"),]$EN.p0, type = "l", lwd=2,
         xlab = xlab, ylab = ylab,
         main = main,
         ylim = c(min(x$EN.p0),max(x$EN.p0)),
         ...)
    lines(x[!is.na(x$INTERIM),]$N, x[!is.na(x$INTERIM),]$EN.p0, type = "o",col="red",lwd=2)
  }
  points(x[x$MIN == "Minimax", "N"], x[x$MIN == "Minimax", "EN.p0"], pch = "M",col="blue")
  points(x[x$OPT == "Optimal", "N"], x[x$OPT == "Optimal", "EN.p0"], pch = "O",col="blue")
  points(x[x$ADMISS == "Admissible", "N"], x[x$ADMISS == "Admissible", "EN.p0"], pch = "A",col="green")
  if (sum(x$EN_opt!="EN.p0_optimal")>0){
    legend("topleft",col=c("black","red"),c("Minimum EN.p0",paste0("Desired Interim at ",unique(x[!is.na(x$INTERIM),]$INTERIM))),lty=1,bty = "n")
  }
}
