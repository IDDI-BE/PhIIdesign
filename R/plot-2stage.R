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
#'                           eps = 0.005, N_min = 0, N_max = 50)
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
  plot(x$N, x$EN.p0, type = "l",
       xlab = xlab, ylab = ylab,
       main = main, ...)
  points(x[x$MIN == "Minimax", "N"], x[x$MIN == "Minimax", "EN.p0"], pch = "M")
  points(x[x$OPT == "Optimal", "N"], x[x$OPT == "Optimal", "EN.p0"], pch = "O")
  points(x[x$ADMISS == "Admissible", "N"], x[x$ADMISS == "Admissible", "EN.p0"], pch = "A")
}

