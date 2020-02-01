#' @export
density.stlpp <- function(X,lbw,tbw,at=c("points","pixels"),dimt=512,...){
  
  if (!inherits(X, "stlpp")) stop("X should an object of class stlpp")
  
  if(missing(at)) at <- "pixels"
  ox <- X$data$x
  oy <- X$data$y
  ot <- X$data$t
  
  L <-  X$domain
  
  n <- npoints(X) # Emerge number of points
  
  stint <- 0 # define the vacant vector to save the densities values per points
  
  if (missing(tbw)) {
    d <- density(ot)
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
  
  pX <- as.stlpp.lpp(X)
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
  if(at=="points") class(out) <- c("numeric","stlppint")
  if(at=="pixels") class(out) <- c("list","stlppint")
  
  attr(out,"tempden") <- d
  attr(out,"netint") <- ldens
  attr(out,"time") <- ot
  attr(out,"bw") <- c("sigma_l"=lbw,"sigma_t"=d$bw)
  
  return(out)
  
}
