#' @export
"[.stlpp" <- function(x, i) {
  d <- as.data.frame(x$data[i,])
  out <- as.stlpp(d$x,d$y,d$t,L=x$domain)
  out$time <- x$time
  return(out)
}


#' @export
"[.stlppint" <- function(x, i){
  
  stopifnot(any(class(i)=="stlpp", class(i)=="numeric"))
  
  if(inherits(i, "stlpp")){
    
    if(!is.null(attr(x,"tgrid"))){
      tgrid <- attr(x,"tgrid")  
    }else{
      tgrid <- attr(x,"tempden")$x
    }
    
    t <- i$data$t
    n <- npoints(i)
    is <- as.lpp.stlpp(i)
    
    id <- findInterval(t,tgrid)
    id[which(id==0)] <- 1
    out <- c()
    for (i in 1:n){
      out[i] <- x[[(id[i])]][is[i]]
    }
    return(out)
  }
  else{
    stlpp <- attr(x,"stlpp")
    return(x[stlpp][i])
  }
}



#' @export
as.data.frame.sumstlpp <- function(x,...){

  r <- rep(x$r,times=attr(x,"nxy"))
  t <- rep(x$t,each=attr(x,"nxy"))
  f <- as.vector(x[[1]])
  ftheo <- as.vector(x[[2]])
  out <- data.frame(r=r,t=t,f=f,ftheo=ftheo,...)
  if(any(names(x)=="Kest"))  colnames(out) <- c("r","t","Kest","Ktheo")
  if(any(names(x)=="Kinhom"))  colnames(out) <- c("r","t","Kinhom","Ktheo")
  if(any(names(x)=="gest"))  colnames(out) <- c("r","t","gest","gtheo")
  if(any(names(x)=="ginhom"))  colnames(out) <- c("r","t","ginhom","gtheo")

  return(out)
}

#' @export
stlpp <- function(X,L,T,...){

    if(missing(L) & !any(class(X)=="lpp")) stop("L is not introduced")
  
  if(!any(class(X)=="lpp")){
    stopifnot(inherits(L, "linnet"))
    Y <- lpp(X,L,...)
    }
  else{
    Y <- X
    L <- domain(X)
    }
  d <- cbind(as.data.frame(Y),t=T)
  
  out <- ppx(data=d[,c(1,2,5)],domain = L,coord.type = c("s","s","t"))
  class(out) <- c("stlpp","ppx")
  out$time <- c(floor(min(T)),ceiling(max(T)))
  return(out)
  
}

#' @export
as.linim.stlppint <- function(X,...){
  if (inherits(X, "stlppint") == FALSE) stop(" X must be from class stlppint")
  if(!is.null(attr(X,"tgrid"))){
    delta <- attr(X,"tgrid")[2]-attr(X,"tgrid")[1]
  }else{
    delta <- attr(X,"tempden")$x[2]-attr(X,"tempden")$x[1]
  }
  
  v <- X[[1]]$v
  v[!is.na(v)] <- 0
  L <- attr(X,"stlpp")$domain
  out <- as.linim(v,L=L,...)
  for (i in 1:length(X)){
    out <- out+X[[i]]
  }
  out <- out*delta
  return(out)
}

#' @export
as.tppint.stlppint <- function(x){
  if (inherits(x, "stlppint") == FALSE) stop(" x must be from class stlppint")
  if(!is.null(attr(x,"tgrid"))){
    delta <- attr(x,"tgrid")
  }else{
    delta <- attr(x,"tempden")$x
  }
  out <- unlist(lapply(x, integral.linim))
  
  class(out) <- c("tppint")
  
  if(!is.null(attr(x,"tgrid"))){
    attr(out,"tgrid") <- attr(x,"tgrid")
  }else{
    attr(out,"tempden") <- attr(x,"tempden")
  }  
  attr(out,"time") <- attr(x,"time")
  
  return(out)
}