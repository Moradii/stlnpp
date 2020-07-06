#' @export
densityVoronoi.tpp <- function(X, f = 1, nrep = 1,
                               at=c("points","pixels"),
                               dimt=128,...){
  
  if(!inherits(X, "tpp")) stop("X should an object of class tpp")
  
  if(missing(at)) at <- "pixels"
  
  n <- npoints(X)
  
  if(f<0 | f>1) stop("f should be between 0 and 1")
  
  Xt <- lpp(X=cbind(X$data$t,rep(0,n)), 
            L=linnet_interval(startp=X$time[1], endp=X$time[2]))
  
  out <- densityVoronoi.lpp(Xt,f=f,nrep=nrep,dimyx=dimt,...)
  
  if(at=="pixels"){
    out1 <- out$v[!is.na(out$v)]
    class(out1) <- c("tppint")
    attr(out1,"time") <- X$data$t
    
    attr(out1,"tgrid") <- out$xcol
    attr(out1,"tpp") <- X
    
    return(out1)
    }
  else{
    Tint <- out$v[!is.na(out$v)]
    ID <- findInterval(X$data$t, out$xcol)
    ID[which(ID==0)] <- 1
    Tintout <- Tint[ID]
    
    class(Tintout) <- c("numeric")
    attr(out,"time") <- X$data$t
    
    return(Tintout)
  }
}


