#' @import graphics
#' @export
plot.stlpp=function(x,xlab=xlab,...){

  if (inherits(x, "stlpp") == FALSE) stop("x should be an object of stlpp")
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar)) 
  
  par(mfrow = c(1, 2), pty = "s")
  plot(as.lpp.stlpp(x), main = "xy-locations on linear network",...)
  
  xx = sort(as.data.frame(x$data[,3])[,1], index.return = TRUE)
  x1 = x$data[xx$ix, ]
  x1=as.data.frame(x1)
  if (missing(xlab)) xlab <- "time"
  plot(x1[, 3], cumsum(x1[, 3]), type = "l", xlab=xlab,
       ylab = "", main = "cumulative number", las = 1
       ,xlim=c(min(x1[,3]),max(x1[,3])),...)
}

#' @export
plot.stlppint <- function(x,style=style,xlab=xlab,xlim=xlim,...){
  
  if (inherits(x, "stlppint") == FALSE) stop(" x must be from class stlppint")
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  t <- attr(x,"time")
  par(mfrow=c(1,2))
  
  if(is.null(attr(x,"tempden")) & is.null(attr(x,"tint"))){
    onL <- as.linim.stlppint(x)
    onT <- as.tppint.stlppint(x)
    
    tgrid <- attr(onT,"tgrid")
    
    if (missing(xlim)) xlim <- range(tgrid)
    if (missing(xlab)) xlab <- "time"
    
    plot(tgrid,onT,ylab="",main="",
         type="l",ylim = c(0,max(onT,table(round(t)))),
         xlab=xlab,xlim = xlim,...)
    points(table(round(t)))
    title(ylab=expression(hat(lambda)[time]), line=2,cex=3,...)
    
    if (missing(style)) {plot(onL,main="",...)}
    else {plot(onL,main="",style=style,...)}
  }
  else if(!is.null(attr(x,"tempden"))){
    d <- attr(x,"tempden")
    int <- length(t)*d$y
    
    if (missing(xlim)) xlim <- range(d$x)
    OK <- d$x>=range(xlim)[1] & d$x<=range(xlim)[2]
    if (missing(xlab)) xlab <- "time"
    plot(d$x[OK],int[OK],
         ylab="",main="",type="l",ylim = c(0,max(int,table(round(t)))),xlab=xlab,xlim = xlim,...)
    points(table(round(t)))
    title(ylab=expression(hat(lambda)[time]), line=2,cex=3,...)
    
    if (missing(style)) {plot(attr(x,"netint"),main="",...)}
    else {plot(attr(x,"netint"),main="",style=style,...)}
  }
  else{
    d <- attr(x,"tint")
    tgrid <- attr(x,"tgrid")
    
    if (missing(xlim)) xlim <- range(tgrid)
    if (missing(xlab)) xlab <- "time"
    
    plot(tgrid,d,ylab="",main="",
         type="l",ylim = c(0,max(d,table(round(t)))),
         xlab=xlab,xlim = xlim,...)
    points(table(round(t)))
    title(ylab=expression(hat(lambda)[time]), line=2,cex=3,...)
    
    if (missing(style)) {plot(attr(x,"lint"),main="",...)}
    else {plot(attr(x,"lint"),main="",style=style,...)}
  }
  
}

#' @export
plot.sumstlpp <- function(x,style=c("level","contour","perspective"),theta=35,phi=10,
                          facets=FALSE,ticktype= "detailed",resfac=5,
                          xlab="r = distance",ylab="t = time",
                          ...){

  if(missing(style)) style="level"
  if (!inherits(x, "sumstlpp")) stop("x is not an object of class sumstlpp")

  if (!requireNamespace("lattice", quietly = TRUE))
    stop("lattice required: install lattice and try again")
  
  if (!requireNamespace("graphics", quietly = TRUE))
    stop("graphics required: install graphics and try again")
  
  if (!requireNamespace("plot3D", quietly = TRUE))
    stop("plot3D required: install plot3D and try again")
  
  if(any(names(x)=="Kest")){

    if (style=="level"){
      plot(lattice::levelplot(x$Kest,row.values=x$r,column.values=x$t,xlab=xlab,ylab=ylab,...))
    }

    if (style=="contour"){
      graphics::contour(x$Kest,x=x$r,y=x$t,xlab=xlab,ylab=ylab,...)
    }

    if (style=="perspective"){
      plot3D::persp3D(x$r,x$t,x$Kest,
              xlab=xlab,ylab=ylab,zlab="",
              main=expression(italic({hat(K)[L]^{ST}}(r,t))),
              zlim= range(c(min(x$Kest,x$Ktheo),max(x$Kest,x$Ktheo))),
              theta=theta,phi=phi,
              facets=facets,ticktype= ticktype,resfac=resfac,
              ...)
    }

  }

  else if (any(names(x)=="gest")){

    if (style=="level"){
      plot(lattice::levelplot(x$gest,row.values=x$r,column.values=x$t,
                              xlab=xlab,ylab=ylab,...))
    }

    if (style=="contour"){
      graphics::contour(x$gest,x=x$r,y=x$t,xlab=xlab,ylab=ylab,...)
    }

    if (style=="perspective"){
      plot3D::persp3D(x$r,x$t,x$gest,
            xlab=xlab,ylab=ylab,zlab="",
            main=expression(italic({hat(g)[L]^{ST}}(r,t))),
            zlim= range(c(min(x$gest,x$gtheo),max(x$gest,x$gtheo))),
            theta=theta,phi=phi,
            facets=facets,ticktype= ticktype,resfac=resfac,
            ...)
   }
    }

  else if (any(names(x)=="Kinhom")){

    if (style=="level"){
      plot(lattice::levelplot(x$Kinhom,row.values=x$r,column.values=x$t,
                              xlab=xlab,ylab=ylab,...))
    }

    if (style=="contour"){
      graphics::contour(x$Kinhom,x=x$r,y=x$t,xlab=xlab,ylab=ylab,...)
    }



    if (style=="perspective"){
      plot3D::persp3D(x$r,x$t,x$Kinhom,
            xlab=xlab,ylab=ylab,zlab="",
            main=expression(italic({hat(K)[LI]^{ST}}(r,t))),
            zlim= range(c(min(x$Kinhom,x$Ktheo),max(x$Kinhom,x$Ktheo))),
            theta=theta,phi=phi,
            facets=facets,ticktype= ticktype,resfac=resfac,
            ...)

    }
    }

  else if (any(names(x)=="ginhom")){

    if (style=="level"){
      plot(lattice::levelplot(x$ginhom,row.values=x$r,column.values=x$t
                              ,xlab=xlab,ylab=ylab,...))
    }

    if (style=="contour"){
      graphics::contour(x$ginhom,x=x$r,y=x$t,xlab=xlab,ylab=ylab,...)
    }

    if (style=="perspective"){
      plot3D::persp3D(x$r,x$t,x$ginhom,
            xlab=xlab,ylab=ylab,zlab="",
            main=expression(italic({hat(g)[LI]^{ST}}(r,t))),
            zlim= range(c(min(x$ginhom,x$gtheo),max(x$ginhom,x$gtheo))),
            theta=theta,phi=phi,
            facets=facets,ticktype= ticktype,resfac=resfac,
            ...)

    }
    }

}
