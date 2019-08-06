#' @export
plot.stlpp=function(X,xlab=xlab,...){

  if (inherits(X, "stlpp") == TRUE) {
    par(mfrow = c(1, 2), pty = "s")
    plot(as.stlpp.lpp(X), main = "xy-locations on linear network",...)
  }

  else (stop("X should be an object of stlpp"))

  xx = sort(as.data.frame(X$data[,3])[,1], index.return = TRUE)
  x = X$data[xx$ix, ]
  x=as.data.frame(x)
  if (missing(xlab)) xlab <- "time"
  plot(x[, 3], cumsum(x[, 3]), type = "l", xlab=xlab,
       ylab = "", main = "cumulative number", las = 1
       ,xlim=c(min(x[,3]),max(x[,3])),...)
}

#' @export
plot.stlppint <- function(X,style=style,xlab=xlab,xlim=xlim,...){
  if (inherits(X, "stlppint") == FALSE) stop(" X must be from class stlppint")
  t <- attr(X,"time")
  par(mfrow=c(1,2))
  d <- attr(X,"tempden")
  int <- length(t)*d$y

  if (missing(xlim)) xlim <- range(d$x)
  OK <- d$x>=range(xlim)[1] & d$x<=range(xlim)[2]
  if (missing(xlab)) xlab <- "time"
  plot(d$x[OK],int[OK],
       ylab="",main="",type="l",ylim = c(0,max(int,table(round(t)))),xlab=xlab,xlim = xlim,...)
  points(table(round(t)))
  title(ylab=expression(hat(lambda)[time]), line=2,cex=3,...)

  if (missing(style)) {plot(attr(X,"netint"),main="",...)}
  else {plot(attr(X,"netint"),main="",style=style,...)}

}

#' @export
plot.sumstlpp <- function(X,style=c("level","contour","perspective"),theta=35,phi=10,
                          facets=FALSE,ticktype= "detailed",resfac=5,
                          xlab="r = distance",ylab="t = time",
                          ...){

  if(missing(style)) style="level"
  if (!inherits(X, "sumstlpp")) stop("X is not an object of class sumstlpp")

  if (!requireNamespace("lattice", quietly = TRUE))
    stop("lattice required: install lattice and try again")
  
  if (!requireNamespace("graphics", quietly = TRUE))
    stop("graphics required: install graphics and try again")
  
  if (!requireNamespace("plot3D", quietly = TRUE))
    stop("plot3D required: install plot3D and try again")
  
  if(any(names(X)=="Kest")){

    if (style=="level"){
      plot(lattice::levelplot(X$Kest,row.values=X$r,column.values=X$t,xlab=xlab,ylab=ylab,...))
    }

    if (style=="contour"){
      graphics::contour(X$Kest,x=X$r,y=X$t,xlab=xlab,ylab=ylab,...)
    }

    if (style=="perspective"){
      plot3D::persp3D(X$r,X$t,X$Kest,
              xlab=xlab,ylab=ylab,zlab="",
              main=expression(italic({hat(K)[L]^{ST}}(r,t))),
              zlim= range(c(min(X$Kest,X$Ktheo),max(X$Kest,X$Ktheo))),
              theta=theta,phi=phi,
              facets=facets,ticktype= ticktype,resfac=resfac,
              ...)
    }

  }

  else if (any(names(X)=="gest")){

    if (style=="level"){
      plot(lattice::levelplot(X$gest,row.values=X$r,column.values=X$t,
                              xlab=xlab,ylab=ylab,...))
    }

    if (style=="contour"){
      graphics::contour(X$gest,x=X$r,y=X$t,xlab=xlab,ylab=ylab,...)
    }

    if (style=="perspective"){
      plot3D::persp3D(X$r,X$t,X$gest,
            xlab=xlab,ylab=ylab,zlab="",
            main=expression(italic({hat(g)[L]^{ST}}(r,t))),
            zlim= range(c(min(X$gest,X$gtheo),max(X$gest,X$gtheo))),
            theta=theta,phi=phi,
            facets=facets,ticktype= ticktype,resfac=resfac,
            ...)
   }
    }

  else if (any(names(X)=="Kinhom")){

    if (style=="level"){
      plot(lattice::levelplot(X$Kinhom,row.values=X$r,column.values=X$t,
                              xlab=xlab,ylab=ylab,...))
    }

    if (style=="contour"){
      graphics::contour(X$Kinhom,x=X$r,y=X$t,xlab=xlab,ylab=ylab,...)
    }



    if (style=="perspective"){
      plot3D::persp3D(X$r,X$t,X$Kinhom,
            xlab=xlab,ylab=ylab,zlab="",
            main=expression(italic({hat(K)[LI]^{ST}}(r,t))),
            zlim= range(c(min(X$Kinhom,X$Ktheo),max(X$Kinhom,X$Ktheo))),
            theta=theta,phi=phi,
            facets=facets,ticktype= ticktype,resfac=resfac,
            ...)

    }
    }

  else if (any(names(X)=="ginhom")){

    if (style=="level"){
      plot(lattice::levelplot(X$ginhom,row.values=X$r,column.values=X$t
                              ,xlab=xlab,ylab=ylab,...))
    }

    if (style=="contour"){
      graphics::contour(X$ginhom,x=X$r,y=X$t,xlab=xlab,ylab=ylab,...)
    }

    if (style=="perspective"){
      plot3D::persp3D(X$r,X$t,X$ginhom,
            xlab=xlab,ylab=ylab,zlab="",
            main=expression(italic({hat(g)[LI]^{ST}}(r,t))),
            zlim= range(c(min(X$ginhom,X$gtheo),max(X$ginhom,X$gtheo))),
            theta=theta,phi=phi,
            facets=facets,ticktype= ticktype,resfac=resfac,
            ...)

    }
    }

}
