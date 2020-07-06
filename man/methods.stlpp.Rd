\name{methods.stlpp}
\alias{methods.stlpp}
\alias{plot.stlpp}
\alias{plot.stlppint}
\alias{plot.sumstlpp}
\alias{print.stlpp}
\alias{print.stlppint}
\alias{print.sumstlpp}
\alias{[.stlpp}
\alias{[.stlppint}
\alias{as.linim.stlppint}
\alias{as.tppint.stlppint}



\title{
Methods for space-time point patterns on a linear network
}
\description{
Methods for space-time point patterns on a linear network.
}
\usage{
\method{plot}{stlpp}(x,xlab = xlab,...)
\method{plot}{stlppint}(x,style=style,xlab=xlab,xlim=xlim,...)
\method{plot}{sumstlpp}(x,style=c("level","contour","perspective"), theta = 35, phi = 10,
facets = FALSE, ticktype = "detailed", resfac = 5,xlab="r = distance",ylab="t = time",...)
\method{print}{stlpp}(x,...)
\method{print}{stlppint}(x,...)
\method{print}{sumstlpp}(x,...)
\method{[}{stlpp}(x,i)
\method{[}{stlppint}(x,i)
\method{as.linim}{stlppint}(X,...)
\method{as.tppint}{stlppint}(x)
}
\arguments{
\item{x,X}{an object of classes \code{\link{stlpp}}, stlppint or sumstlpp
}
\item{style}{
style of plot
}
\item{theta,phi}{
see \link[plot3D]{persp3D}
}
\item{facets,ticktype}{
see \link[plot3D]{persp3D}
}
\item{resfac}{
see \link[plot3D]{persp3D}
}
\item{xlab,ylab}{
the x,y label of the plot
}
\item{xlim}{
giving the x limits for the plot
}
\item{i}{
numeric, logical, or an object of class \code{\link{stlpp}}
}
\item{...}{
either ignore for \code{\link{as.linim}}, or
gracphical arguments passed to \code{\link{plot}}/\code{\link{print}}
}
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}