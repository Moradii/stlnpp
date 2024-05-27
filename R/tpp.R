#' Create a temporal point pattern
#'
#' Create an object of class \code{\link{tpp}} representing a one-dimensional point pattern.
#'
#' @usage tpp(X,a,b)
#'
#' @param X an object of class \code{\link{numeric}}, \code{\link{integer}} or \code{\link{vector}}
#' @param a lower band of the time domain. if not given by the user, it will be the minimum of X
#' @param b upper bound of the time domain. if not given by the user, it will be the maximum of X
#' 
#' @seealso \code{\link{stlpp}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com>
#' 
#' @returns 
#' An object of class tpp.
#' 
#' @details 
#' Create a one-dimensional point pattern.
#'
#' 
#' @examples  
#' tpp(runif(10))
#'
#' @import spatstat
#' @import spatstat.geom
#' @import spatstat.linnet
#' @importFrom spatstat.univar dkernel
#' @import stats
#' @export
tpp <- function(X,a,b){
  
  stopifnot(inherits(X,"numeric") | inherits(X,"integer") | inherits(X,"vector"))
  out <- ppx(data=X,coord.type = c("t"))
  names(out$data) <- "t"
  class(out) <- c("tpp","ppx")
  if(missing(a)) a <- floor(min(X))
  if(missing(b)) b <- ceiling(max(X))
  out$time <- c(a,b)
  return(out)
  
}