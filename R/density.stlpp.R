#' @import spatstat
#' @import stats
#' @export
density.stlpp <- function(x,lbw,tbw,at=c("points","pixels"),dimt=512,...){
  
  if (!inherits(x, "stlpp")) stop("x should an object of class stlpp")
  
  if(missing(at)) at <- "pixels"
  ox <- x$data$x
  oy <- x$data$y
  ot <- x$data$t
  
  L <-  x$domain
  
  n <- npoints(x) # Emerge number of points
  
  stint <- 0 # define the vacant vector to save the densities values per points
  
  if (missing(tbw)) {
    d <- density(ot,n=dimt)
  }
  else{
    d <- density(ot,bw=tbw,n=dimt)
  }
  
  if(at=="points"){
    Tint <- d$y[findInterval(ot, d$x)] * n
  }
  else{
    Tint <- d$y * n
  }
  ############################################## space intensity
  
  pX <- as.lpp.stlpp(x)
  if (missing(lbw)) lbw <- bw.scott.iso(pX)
  
  if(at=="points"){
    ldens <- density.lpp(pX,sigma = lbw,distance="euclidean",...)
    Sint <- density.lpp(pX,sigma = lbw,at="points",distance="euclidean",...)
    stint <-  Sint*Tint/npoints(pX)
    out <- stint
  }
  
  if(at=="pixels"){
    
    ldens <- density.lpp(pX,sigma = lbw,distance="euclidean",...)
    out <- lapply(X=1:length(Tint), function(i){
      ldens*Tint[i]/npoints(pX)
    })
  }
  
  
  names(lbw) <- NULL
  if(at=="points") class(out) <- c("numeric")
  if(at=="pixels") class(out) <- c("list","stlppint")
  
  attr(out,"tempden") <- d
  attr(out,"netint") <- ldens
  attr(out,"time") <- ot
  attr(out,"bw") <- c("sigma_l"=lbw,"sigma_t"=d$bw)
  attr(out,"stlpp") <- x
  
  
  return(out)
  
}
