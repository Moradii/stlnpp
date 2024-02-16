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
print.tppint <- function(x,...){
  print(as.vector(x),...)
}

#' @export
plot.tppint <- function(x,xlab=xlab,xlim=xlim,line=2.5,main=NULL,...){
  if (inherits(x, "tppint") == FALSE) stop(" x must be from class tppint")
  t <- attr(x,"time")
  
  if(!is.null(attr(x,"tempden"))){
    d <- attr(x,"tempden")
    int <- length(t)*d$y
    tgrid <- d$x
  }else{
    tgrid <- attr(x,"tgrid")
    int <- x 
  }
  
  
  if (missing(xlim)) xlim <- range(tgrid)
  
  OK <- tgrid>=range(xlim)[1] & tgrid<=range(xlim)[2]
  
  if (missing(xlab)) xlab <- "time"
  
  plot(tgrid[OK],as.numeric(int)[OK],
       ylab="",main=main,type="l",ylim = c(0,max(int,table(round(t)))),xlab=xlab,xlim = xlim,...)
  points(table(round(t)))
  title(ylab=expression(hat(lambda)[time]), line=line,cex=3,...)
}


#' @export
"[.tpp" <- function(x, i) {
  stopifnot(any(class(i)=="numeric", class(i)=="logical"))
  
  d <- as.data.frame(x$data[i,])
  out <- tpp(d$t)
  out$time <- x$time
  return(out)
}


#' @export
"[.tppint" <- function(x, i){
  
  stopifnot(any(class(i)=="tpp", class(i)=="numeric", class(i)=="logical"))
  
  if(inherits(i, "tpp")){
    
    if(!is.null(attr(x,"tgrid"))){
      tgrid <- attr(x,"tgrid")  
    }else{
      tgrid <- attr(x,"tempden")$x
    }
    
    t <- i$data$t
    n <- npoints(i)
    # is <- as.lpp.stlpp(i)
    
    id <- findInterval(t,tgrid)
    id[which(id==0)] <- 1
    out <- c()
    for (j in 1:n){
      out[j] <- as.numeric(x)[id[j]]
    }
    return(out)
  }
  else{
    tp <- attr(x,"tpp")
    return(x[tp][as.numeric(i)])
  }
}