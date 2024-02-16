#' Methods for spatio-temporal point patterns on a linear network
#'
#' This function projects an object of class \code{\link{stlpp}} into a linear network.
#'
#' @usage \method{as.lpp}{stlpp}(x,...)
#'
#' @param x an object of class \code{\link{stlpp}}
#' @param ... arguments passed to \code{\link{as.lpp}}
#' 
#' @seealso \code{\link{as.stlpp}}, \code{\link{lpp}}, \code{\link{as.lpp}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com>
#' 
#' @returns An object of class \code{\link{lpp}}.
#' 
#' @details This function projects the spatio-temporal point pattern x on the linear network L into L, giving its corresponding spatial point pattern.
#' 
#' @examples  
#' data(easynet)
#' x <- runifpointOnLines(40, easynet)
#' t1 <- sample(1:10,40,replace=TRUE)
#' Y <- as.stlpp(x,t=t1,L=easynet)
#' as.lpp.stlpp(Y)
#' 
#' @export
as.lpp.stlpp <- function(x,...){
  if(!any(class(x)=='stlpp')) stop('class(x) must be stlpp')
  Y <- as.lpp(x=x$data$x,y=x$data$y,L=x$domain,...)
  return(Y)
}
