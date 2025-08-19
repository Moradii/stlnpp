#' Methods for spatio-temporal point patterns on a linear network
#'
#' This function projects an object of class \code{\link{stlpp}} into a linear network.
#'
#' @usage \method{as.lpp}{stlpp}(x,...)
#'
#' @param x an object of class \code{\link{stlpp}}
#' @param ... arguments passed to \code{\link[spatstat.linnet]{as.lpp}}
#' 
#' @seealso \code{\link{as.stlpp}}, \code{\link[spatstat.linnet]{lpp}}, \code{\link[spatstat.linnet]{as.lpp}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com>
#' 
#' @returns An object of class \code{\link[spatstat.linnet]{lpp}}.
#' 
#' @details This function projects the spatio-temporal point pattern x into its corresponding linear network L which is the corresponding spatial point pattern of x
#' 
#' @examples  
#' data(easynet)
#' x <- runifpointOnLines(40, easynet)
#' t1 <- sample(1:10, 40, replace = TRUE)
#' Y <- as.stlpp(x, t = t1, L = easynet)
#' as.lpp.stlpp(Y)
#' 
#' @export
as.lpp.stlpp <- function(x,...){
  if(!any(class(x)=='stlpp')) stop('class(x) must be stlpp')
  Y <- spatstat.linnet::as.lpp(x=x$data$x,y=x$data$y,L=x$domain,...)
  return(Y)
}
