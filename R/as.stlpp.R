#' Convert data to a spatio-temporal point pattern on a linear network
#'
#' This function converts data to a spatio-temporal point pattern on a linear network.
#'
#' @usage as.stlpp(x,y,t,L)
#'
#' @param x,y,t vectors of Cartesian coordinates and time occurrence. Alternatively, x can be of classes \code{\link{data.frame}}, \code{\link{ppp}} and \code{\link{lpp}}
#' @param L linear network (object of class \code{\link{linnet}})
#' 
#' @seealso \code{\link{stlpp}}
#' 
#' @author Mehdi Moradi <m2.moradi@yahoo.com>
#' 
#' @returns A spatio-temporal point pattern on a linear network. An object of class  \code{\link{stlpp}}.
#' 
#' @details This function converts data to an object of class stlpp.
#'  Data can be of formats:
#'  \itemize{
#'    \item x is of class class \code{\link{data.frame}} with three columns. Then columns are considered as Cartesian coordinates (i.e. x,y,t) and they will be converted to a spatio-temporal point pattern on the linear network L.
#'    
#'    \item x is a planar point pattern (class \code{\link{ppp}}). Then x will be converted to a spatio-temporal point pattern on the linear network L and with coresponding time vector t.
#'    
#'    \item x is a linear point pattern (class \code{\link{lpp}}). Then x will be converted to a spatio-temporal point pattern on the linear network L and with coresponding time vector t.
#'    
#'    \item x,y,t are vectors of same length where x,y are living on the corresponding network L.
#'  }
#' 
#' @examples  
#' data(easynet)
#' x <- runifpointOnLines(40, easynet)
#' t1 <- sample(1:10,40,replace=TRUE)
#' Y <- as.stlpp(x,t=t1,L=easynet)
#'
#' Z <- as.lpp.stlpp(Y)
#' t2 <- sample(1:10,40,replace=TRUE)
#' W <- as.stlpp(Z,t=t2)
#' 
#' @export
as.stlpp <- function(x=NULL,y=NULL,t=NULL,L=NULL){

  if(any(class(x)=="ppp")){
    data <- data.frame(as.data.frame(x),t)
    #colnames(data) <- c("x","y","t")

    if(missing(L)) stop("L is not introduced")
    out <- ppx(data=data,domain = L,coord.type = c("s","s","t"))
    class(out) <- c("stlpp","ppx")
    out$time <- range(t)
    return(out)
  }

  if(any(class(x)=="lpp")){
    data <- data.frame(as.data.frame(x)[,1:2],t)
    #colnames(data) <- c("x","y","t")
    out <- ppx(data=data,domain = domain(x),coord.type = c("s","s","t"))
    class(out) <- c("stlpp","ppx")
    out$time <- range(t)
    return(out)
  }

  if(any(class(x)=="data.frame")){
    if(missing(L)) stop("L is not introduced")
    data <- x
    colnames(data) <- c("x","y","t")
    out <- ppx(data=data,domain = L,coord.type = c("s","s","t"))
    class(out) <- c("stlpp","ppx")
    out$time <- range(data$t)
    return(out)
  }

  if(missing(L)) stop("L is not introduced")
  if(!(length(x)==length(y) & length(y)==length(t))) stop("x,y,t are not of same length.")
  pointonl <- lpp(list(x=x,y=y),L=L)
  data <- data.frame(pointonl$data$x,pointonl$data$y,t)
  colnames(data) <- c("x","y","t")
  out <- ppx(data=data,domain = L,coord.type = c("s","s","t"))
  class(out) <- c("stlpp","ppx")
  out$time <- range(t)
  return(out)
}

