#' Create spatio-temporal point pattern on linear network
#'
#' Create an object of class \code{\link{stlpp}} representing a spatio-temporal point pattern on a linear network.
#'
#' @usage stlpp(X, L, T, ...)
#'
#' @param X Locations of the points. a matrix or data frame of coordinates, or a point pattern object (of class "ppp") or other data acceptable to \code{\link{as.ppp}} or \code{\link{lpp}}
#' @param L linear network (object of class \code{\link{linnet}}) on which the points lie
#' @param T time occurrence of the points
#' @param ... ignored
#' 
#' 
#' @seealso \code{\link{as.stlpp}}, \code{\link{lpp}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com>
#' 
#' @returns An object of class  \code{\link{stlpp}}.
#' 
#' @details This function creates an object of class \code{\link{stlpp}}. For details about X see \code{\link{lpp}}. \code{T} represents the time occurrences of data points.
#' 
#' @examples  
#' data(easynet)
#' X <- rpoislpp(1,easynet)
#' t <- runif(npoints(X))
#' stlpp(X,T=t,L=easynet)
#' 
#' @export
stlpp <- function(X,L,T,...){
  
  if(missing(L) & !any(class(X)=="lpp")) stop("L is not introduced")
  
  if(!any(class(X)=="lpp")){
    stopifnot(inherits(L, "linnet"))
    Y <- lpp(X,L,...)
  }
  else{
    Y <- X
    L <- domain(X)
  }
  d <- cbind(as.data.frame(Y),t=T)
  
  out <- ppx(data=d[,c(1,2,5)],domain = L,coord.type = c("s","s","t"))
  class(out) <- c("stlpp","ppx")
  out$time <- c(floor(min(T)),ceiling(max(T)))
  return(out)
  
}