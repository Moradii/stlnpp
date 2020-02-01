#' @export
tpp <- function(X){
  
  stopifnot(inherits(X,"numeric") | inherits(X,"integer") | inherits(X,"vector"))
  out <- ppx(data=X,coord.type = c("t"))
  names(out$data) <- "t"
  class(out) <- c("tpp","ppx")
  out$time <- round(range(X),4)
  return(out)
  
}

#' @export
as.stlpp.tpp <- function(X){
  if(!any(class(X)=="stlpp")) stop("class(X) must be stlpp")
  out <- ppx(data=X$data$t,coord.type = c("t"))
  names(out$data) <- "t"
  class(out) <- c("tpp","ppx")
  out$time <- X$time
  return(out)
}

#' @export
print.tpp <- function(x)
{
  if(!any(class(x)=="tpp")) stop("class(X) must be tpp")
  cat("Temporal point pattern \n");
  if(npoints(x)>1){cat(paste0(npoints(x)," ", "points"),"\n")}
  else{cat(paste0(npoints(x)," ", "point"),"\n")};
  cat(paste0("Time period: [",range(x$time)[1],", ", range(x$time)[2],"]"),"\n")
}


#' @export
plot.tpp <- function(X,xlab="time",ylab="",main = "cumulative number",...){
  if(!any(class(X)=="tpp")) stop("class(X) must be tpp")
  xx <-  sort(X$data$t, index.return = TRUE)
  x  <-  X$data$t[xx$ix]
  plot(x, cumsum(x), type = "l", las = 1,xlab=xlab,ylab=ylab,main=main,
       xlim=c(min(x),max(x)),...)
}

#' @export
density.tpp <- function(X,tbw,at=c("points","pixels"),...){
  
  if (!inherits(X, "tpp")) stop("X should an object of class tpp")
  
  n <- npoints(X) # Emerge number of points
  
  if (missing(tbw)) {d <- density(X$data$t,...)}
  else {
    d <- density(X,bw=tbw,...)
  }
  
  if(missing(at)) at <- "pixels"
  
  if(at=="points"){
    Tint <- d$y[findInterval(X$data$t, d$x)] * npoints(X)
  }
  if(at=="pixels"){
    Tint <- d$y * npoints(X)
  }

  out <- Tint
  
  attr(out,"tempden") <- d
  attr(out,"bw") <- d$bw
  attr(out,"time") <- X$data$t
  class(out) <- c("tppint")
  return(out)
}

#' @export
print.tppint <- function(X){
  print(as.vector(X))
}

#' @export
plot.tppint <- function(X,style=style,xlab=xlab,xlim=xlim,line=2.5,...){
  if (inherits(X, "tppint") == FALSE) stop(" X must be from class tppint")
  t <- attr(X,"time")
  d <- attr(X,"tempden")
  int <- length(t)*d$y
  
  if (missing(xlim)) xlim <- range(d$x)
  OK <- d$x>=range(xlim)[1] & d$x<=range(xlim)[2]
  if (missing(xlab)) xlab <- "time"
  plot(d$x[OK],int[OK],
       ylab="",main="",type="l",ylim = c(0,max(int,table(round(t)))),xlab=xlab,xlim = xlim,...)
  points(table(round(t)))
  title(ylab=expression(hat(lambda)[time]), line=line,cex=3,...)
}
