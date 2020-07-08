#' @import spatstat
#' @import stats
#' @export
densityVoronoi.stlpp <- function(X, f = 1, nrep = 1,
                                 separable=FALSE,at=c("points","pixels"),
                                 dimt=128,...){
  
  if(!inherits(X, "stlpp")) stop("X should an object of class stlpp")
  
  if(missing(at)) at <- "pixels"

  n <- npoints(X)
  
  if(f<0 | f>1) stop("f should be between 0 and 1")
  
  if(f==1) separable <- TRUE
  
  Xt <- lpp(X=cbind(X$data$t,rep(0,n)), 
            L=linnet_interval(startp=X$time[1], endp=X$time[2]))
  
  Xs <- as.lpp.stlpp(X)
  
if(separable){
    
    IntEstL <- densityVoronoi(X=Xs,f=f,nrep=nrep,...)
    IntEstT <- densityVoronoi(X=Xt,f=f,nrep=nrep,dimyx=dimt,...)
    
    tgrid <- IntEstT$xcol
    
    IntEstTv <- IntEstT$v[!is.na(IntEstT$v)]/n
    
    
    out <- lapply(X=(1:length(IntEstTv)), FUN=function(j){IntEstL*IntEstTv[j]}) 
    

}else{
    
Y <- rthin(X,P=f,nsim = nrep)
  
out.nonsep <- lapply(X=1:nrep, function(i){
                     densityVoronoi(Y[[i]],f=1,nrep = 1, separable = TRUE,dimt=dimt,...)
                   })
tgrid <- attr(out.nonsep[[1]],"tgrid")
  
out <- list()
k <- length(out.nonsep[[1]])
for (i in 1:k){
  
    out[[i]] <- out.nonsep[[1]][[1]]
    out[[i]]$v[!is.na(out[[i]]$v)] <- 0
   
for (j in 1:length(out.nonsep)){
      out[[i]] <- out[[i]] + out.nonsep[[j]][[i]]
                                     
                               }
    out[[i]] <- out[[i]]/(f*nrep)
                                     }
    }
  
if(at=="points"){
    t <- X$data$t
    id <- findInterval(t,tgrid)
    out1 <- c()
    for (i in 1:n){
      out1[i] <- out[[id[i]]][Xs[i]]
                  }
    out <- out1
                }

  if(at=="points") class(out) <- c("numeric","stlppint")
  if(at=="pixels") class(out) <- c("list","stlppint")
  
  if(separable){
    attr(out,"lint") <- IntEstL
    attr(out,"tint") <- IntEstT$v[!is.na(IntEstT$v)]
  }
  attr(out,"tgrid") <- tgrid
  attr(out,"time") <- X$data$t
  attr(out,"stlpp") <- X
return(out)
}




linnet_interval <- function(startp=0, endp=1,...){
  Wt <- c(startp, endp)
  vertices <- ppp(x=c(Wt[1],Wt[2]), y=c(0,0), 
                  window=owin(Wt,c(-0.5,0.5)))
  m <- matrix(data=c(FALSE,TRUE,TRUE,FALSE), nrow=2, ncol=2, 
              byrow=TRUE)
  out <- linnet(vertices=vertices, m=m,...)
  return(out)
}
