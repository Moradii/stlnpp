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
\alias{as.data.frame.sumstlpp}



\title{
Methods for space-time point patterns on a linear network
}
\description{
Methods for space-time point patterns on a linear network.
}
\usage{
\method{plot}{stlpp}(x,
  onlypoints=FALSE,
  onlytime=FALSE,
  args.xy = list(),
  args.time = list(),
  ...)
\method{plot}{stlppint}(x, style = NULL, xlab = NULL, ylab = NULL, main =
                   NULL, xlim = NULL, ylim = NULL, bar = TRUE,
                   args.spatial = list(), args.temporal = list(),
                   all.images = FALSE, per.page = 9, ...)
\method{plot}{sumstlpp}(x,
                          style = c( "perspective",  "contour","level"),
                          theta = 35, phi = 10, facets = FALSE,
                          ticktype = "detailed", resfac = 5,
                          xlab = NULL, ylab = NULL, zlab = NULL, main = NULL,
                          args.level = list(),
                          args.contour = list(),
                          args.perspective = list(),...)
\method{print}{stlpp}(x,...)
\method{print}{stlppint}(x,...)
\method{print}{sumstlpp}(x,...)
\method{[}{stlpp}(x,i)
\method{[}{stlppint}(x,i)
\method{as.linim}{stlppint}(X,...)
\method{as.tppint}{stlppint}(x)
\method{as.data.frame}{sumstlpp}(x,...)
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
\item{xlab,ylab,zlab}{
the x,y,z label of the plot
}
\item{main}{
an overall title for the plot
}
\item{args.spatial}{
A list of additional arguments passed to spatial plots
}
\item{args.temporal}{
A list of additional arguments passed to temporal plots
}
\item{all.images}{
Logical, if \code{TRUE} plots all spatial images over time, where each image corresponds to the estimated spatio-temporal intensities
}
\item{per.page}{
Number of spatial images to show per plotting page when \code{all.images = TRUE}
}
\item{args.level}{
A named list of additional arguments to pass to \code{\link[lattice]{levelplot}}
}
\item{args.contour}{
A named list of additional arguments to pass to \code{\link[graphics]{contour}}
}
\item{args.perspective}{
A named list of additional arguments to pass to \code{\link[plot3D]{persp3D}}
}
\item{xlim}{
giving the x limits for the plot
}
\item{ylim}{
giving the y limits for the plot
}
\item{bar}{
if TRUE, bar plot of rounded time occurances will be added to the density plot
}
\item{i}{
numeric, logical, or an object of class \code{\link{stlpp}}
}
\item{onlypoints, onlytime}{
if only points or times are to be plotted
}
\item{
args.xy, args.time
}{
arguments passed to \code{\link[spatstat.linnet]{plot.lpp}} and \code{\link[base]{plot}}
}
\item{...}{
either ignore for \code{\link[spatstat.linnet]{as.linim}}, or
gracphical arguments passed to \code{\link{plot}}/\code{\link{print}}
}
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}