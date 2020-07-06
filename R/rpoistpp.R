#' @export
rpoistpp <- function(lambda,a=NULL,b=NULL,check=FALSE,lmax=NULL,nsim=1){
  
  
  if (!is.numeric(lambda) & !is.function(lambda) & !class(lambda)=="tppint")
    stop(" lambda should be a number or a function")
  
  if(nsim > 1) {
    out <- list()
    for (i in 1:nsim) {
      out[[i]] <- rpoistpp(lambda,a=a,b=b,check=check,lmax=lmax,nsim=1)
    }
    return(out)
  }
  
  if (class(lambda)=="numeric"){
    
    if (a >= b)  stop("lower bound must be smaller than upper bound")
    
    n <- rpois(1,lambda*(b-a))
    t <- runif(n,a,b)
    out <- tpp(t,a=a,b=b)
  }
  else if(is.function(lambda)){
    
    if (a >= b)  stop("lower bound must be smaller than upper bound")
    
    if(is.null(lmax)){
      tgrid <- seq(a,b,length.out = 128)
      lmax <- max(lambda(tgrid))
    }
    
    mean <- lmax*(b-a)
    
    n <- rpois(1,mean)
    
    t <- runif(n,a,b)
    
    prob <- lambda(t)/lmax
    
    if(check) {
      if(any(prob < 0))
        warning("Negative values of lambda obtained")
      if(any(prob > 1))
        warning("lmax is not an upper bound for lambda")
    }
    
    u <- runif(length(t))
    
    retain <-  (u <= prob)
    
    t <- t[retain]
    out <- tpp(t,a=a,b=b)
    
  }else if(any(class(lambda)=="tppint")){
    
    lmax <- max(lambda)
       
       a <- attr(lambda,"tpp")$time[1]
       
       b <- attr(lambda,"tpp")$time[2]
       
    mean <- lmax*(b-a)
       
       n <- rpois(1,mean)
       
       t <- runif(n,a,b)
       
    prob <- lambda[tpp(t)]/lmax
       
       if(check) {
         if(any(prob < 0))
           warning("Negative values of lambda obtained")
         if(any(prob > 1))
           warning("lmax is not an upper bound for lambda")
       }
       
       u <- runif(length(t))
       
       retain <-  (u <= prob)
       
       t <- t[retain]
       out <- tpp(t,a=a,b=b)
  }
  
  class(out) <- c("tpp","ppx")
  
  out$time <- c(a,b)
  
  return(out)
}