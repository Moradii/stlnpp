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
print.tpp <- function(x,...)
{
  if(!any(class(x)=="tpp")) stop("class(X) must be tpp")
  cat("Temporal point pattern \n");
  if(npoints(x)>1){cat(paste0(npoints(x)," ", "points"),"\n")}
  else{cat(paste0(npoints(x)," ", "point"),"\n")};
  cat(paste0("Time period: [",range(x$time)[1],", ", range(x$time)[2],"]"),"\n")
}


#' @export
plot.tpp <- function(x,xlab="time",ylab="",main = "cumulative number",...){
  if(!any(class(x)=="tpp")) stop("class(x) must be tpp")
  xx <-  sort(x$data$t, index.return = TRUE)
  x1  <-  x$data$t[xx$ix]
  plot(x1, cumsum(x1), type = "l", las = 1,xlab=xlab,ylab=ylab,main=main,
       xlim=c(min(x1),max(x1)),...)
}

#' @export
density.tpp <- function(x,tbw,at=c("points","pixels"),...){
  
  if (!inherits(x, "tpp")) stop("x should an object of class tpp")
  
  n <- npoints(x) # Emerge number of points
  
  if (missing(tbw)) {d <- density(x$data$t,...)}
  else {
    d <- density(x,bw=tbw,...)
  }
  
  if(missing(at)) at <- "pixels"
  
  if(at=="points"){
    Tint <- d$y[findInterval(x$data$t, d$x)] * npoints(x)
  }
  if(at=="pixels"){
    Tint <- d$y * npoints(x)
  }

  out <- Tint
  
  attr(out,"tempden") <- d
  attr(out,"bw") <- d$bw
  attr(out,"time") <- x$data$t
  class(out) <- c("tppint")
  return(out)
}

#' @export
print.tppint <- function(x,...){
  print(as.vector(x),...)
}

#' @export
plot.tppint <- function(x,xlab=xlab,xlim=xlim,line=2.5,...){
  if (inherits(x, "tppint") == FALSE) stop(" x must be from class tppint")
  t <- attr(x,"time")
  d <- attr(x,"tempden")
  int <- length(t)*d$y
  
  if (missing(xlim)) xlim <- range(d$x)
  OK <- d$x>=range(xlim)[1] & d$x<=range(xlim)[2]
  if (missing(xlab)) xlab <- "time"
  plot(d$x[OK],int[OK],
       ylab="",main="",type="l",ylim = c(0,max(int,table(round(t)))),xlab=xlab,xlim = xlim,...)
  points(table(round(t)))
  title(ylab=expression(hat(lambda)[time]), line=line,cex=3,...)
}
