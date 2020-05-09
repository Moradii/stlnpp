#' @export
rthin.stlpp <- function(X,P=P,nsim=1){

  if(!any(class(X)=="stlpp")) stop("class(X) must be stlpp")
  if (P<0 | P>1) stop("P should be a numeric value between 0 and 1")

  if (nsim==1){
    u <- runif(npoints(X),0,1)
    OK <- (u<P)
    return(sub.stlpp(X,OK))

  }else {
    out <- list()
    for (i in 1:nsim) {
      u <- runif(npoints(X),0,1)
      OK <- (u<P)
      out[[i]] <- sub.stlpp(X,OK)
    }
    class(out) <- "list"
    return(out)
  }
}

sub.stlpp <- function(x,i){
  d <- as.data.frame(x$data[i,])
  as.stlpp(d$x,d$y,d$t,L=x$domain)
}