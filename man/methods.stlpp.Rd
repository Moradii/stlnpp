\name{methods.stlpp}
\alias{methods.stlpp}
\alias{plot.stlpp}
\alias{plot.stlppint}
\alias{plot.sumstlpp}
\alias{print.stlpp}
\alias{print.stlppint}
\alias{print.sumstlpp}


\title{
Methods for space-time point patterns on a linear network
}
\description{
Methods for space-time point patterns on a linear network.
}
\usage{
\method{plot}{stlpp}(X,xlab = xlab,...)
\method{plot}{stlppint}(X,style=style,xlab=xlab,xlim=xlim,...)
\method{plot}{sumstlpp}(X,style=c("level","contour","perspective"), theta = 35, phi = 10,
facets = FALSE, ticktype = "detailed", resfac = 5,xlab="r = distance",ylab="t = time",...)
\method{print}{stlpp}(X)
\method{print}{stlpp}(X)
\method{print}{stlpp}(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{an object of classes \code{\link{stlpp}}, stlppint or sumstlpp
}
\item{style}{style of plot.}
\item{theta,phi}{see \link[plot3D]{persp3D}.}
\item{facets,ticktype}{see \link[plot3D]{persp3D}.}
\item{resfac}{see \link[plot3D]{persp3D}.}
\item{xlab,ylab}{the x,y label of the plot.}
\item{xlim}{giving the x limits for the plot.}
\item{...}{ gracphical arguments passed to \code{\link{plot}}.}
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}