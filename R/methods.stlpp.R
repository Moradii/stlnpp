#' @export
"[.stlpp" <- function(X, i, j, drop=FALSE) {
  d <- as.data.frame(X$data[i,])
  out <- as.stlpp(d$x,d$y,d$t,L=X$domain)
  out$time <- X$time
  return(out)
}

#' @export
as.data.frame.sumstlpp <- function(X){

  r <- rep(X$r,times=attr(X,"nxy"))
  t <- rep(X$t,each=attr(X,"nxy"))
  f <- as.vector(X[[1]])
  ftheo <- as.vector(X[[2]])
  out <- data.frame(r=r,t=t,f=f,ftheo=ftheo)
  if(any(names(X)=="Kest"))  colnames(out) <- c("r","t","Kest","Ktheo")
  if(any(names(X)=="Kinhom"))  colnames(out) <- c("r","t","Kinhom","Ktheo")
  if(any(names(X)=="gest"))  colnames(out) <- c("r","t","gest","gtheo")
  if(any(names(X)=="ginhom"))  colnames(out) <- c("r","t","ginhom","gtheo")

  return(out)
}

#' @export
stlpp <- function(X,L,T,...){
  stopifnot(inherits(L, "linnet"))
  if(missing(L)) stop("L is not introduced")
  
  Y <- lpp(X,L,...)
  d <- cbind(as.data.frame(Y),T)
  
  out <- ppx(data=d[,c(1,2,5)],domain = L,coord.type = c("s","s","t"))
  class(out) <- c("stlpp","ppx")
  out$time <- round(range(T),4)
  return(out)
  
}

