#' Convert data to a one-dimensional point pattern
#'
#' This function converts an object of class \code{\link{stlpp}} to  class \code{\link{tpp}}.
#'
#' @usage as.tpp.stlpp(X)
#'
#' @param X an object of class \code{\link{stlpp}}
#' 
#' @seealso \code{\link{as.stlpp}}, \code{\link{lpp}}, \code{\link{as.lpp}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com>
#' 
#' @returns An object of class tpp.
#' 
#' @details This function projects the spatio-temporal point pattern X into its corresponding time domain T.
#' 
#' @examples  
#' X <- rpoistlpp(10,1,2,easynet)
#' as.tpp.stlpp(X)
#' 
#' @export
as.tpp.stlpp <- function(X){
  if(!any(class(X)=="stlpp")) stop("X must be of class stlpp")
  out <- tpp(X$data$t)
  # out <- ppx(data=X$data$t,coord.type = c("t"))
  # names(out$data) <- "t"
  # class(out) <- c("tpp","ppx")
  # out$time <- X$time
  return(out)
}