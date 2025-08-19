#' Convert data to a one-dimensional point pattern
#'
#' This function converts an object of class \code{\link{stlpp}} to  class \code{\link{tpp}}.
#'
#' @usage \method{as.tpp}{stlpp}(X)
#'
#' @param X an object of class \code{\link{stlpp}}
#' 
#' @seealso \code{\link{as.stlpp}}, \code{\link[spatstat.linnet]{lpp}}, \code{\link[spatstat.linnet]{as.lpp}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com>
#' 
#' @returns An object of class  \code{\link{tpp}}.
#' 
#' @details This function projects the spatio-temporal point pattern X into its corresponding time domain T.
#' 
#' @examples  
#' X <- rpoistlpp(10,1,2,easynet)
#' as.tpp(X)
#' 
#' @export
#' @rdname as.tpp
as.tpp.stlpp <- function(X){
  if(!any(class(X)=="stlpp")) stop("X must be of class stlpp")
  out <- tpp(X$data$t, a = X$time[1], b = X$time[2])
  # out <- ppx(data=X$data$t,coord.type = c("t"))
  # names(out$data) <- "t"
  # class(out) <- c("tpp","ppx")
  # out$time <- X$time
  return(out)
}



#' @export
#' @rdname as.tpp
as.tpp <- function(X) {
  UseMethod("as.tpp")
}