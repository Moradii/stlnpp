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
  data <- data.frame(x,y,t)
  colnames(data) <- c("x","y","t")
  out <- ppx(data=data,domain = L,coord.type = c("s","s","t"))
  class(out) <- c("stlpp","ppx")
  out$time <- range(t)
  return(out)
}

