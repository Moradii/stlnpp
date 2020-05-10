\name{densityVoronoi.tpp}
\alias{densityVoronoi.tpp}
\title{
Intensity estimate of temporal point patterns using Voronoi-Dirichlet tessellation
}
\description{
adaptive intensity estimation for temporal point patterns using Voronoi-Dirichlet tessellation.
}
\usage{
\method{densityVoronoi}{tpp}(X, f = 1, nrep = 1,at=c("points","pixels"),
                               dimt=128,...)
}
\arguments{
\item{X}{
an object of class \code{\link{tpp}}
}
\item{f}{
fraction (between 0 and 1 inclusive) of the data points that will be used to build a tessellation for the intensity estimate
}
\item{nrep}{
number of independent repetitions of the randomised procedure
}
\item{at}{
string specifying whether to return the intensity values at a grid of pixel locations and time (at="pixels") or only at the points of X (at="points"). default is to estimate the intensity at pixels
}
\item{dimt}{
the number of equally spaced points at which the temporal density is to be estimated. see \link[stats]{density}
}
\item{...}{
arguments passed to \code{\link{densityVoronoi.lpp}}
}
}
\details{
This function computes intensity estimates for temporal point patterns using Voronoi-Dirichlet tessellation.
}
\value{
If \code{at="points"}: a vector of intensity values at the data points of X.

If \code{at="pixels"}: a vector of intensity values over a grid.
}
\references{
Mateu, J., Moradi, M., & Cronie, O. (2019). Spatio-temporal point patterns on linear networks: Pseudo-separable intensity estimation. Spatial Statistics, 100400.
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>,  and Ottmar Cronie.
}
\seealso{
\code{\link{densityVoronoi.lpp}}, \code{\link{densityVoronoi.stlpp}}
}
\examples{
X <- rpoistlpp(0.2,a=0,b=5,L=easynet)
Y <- as.tpp.stlpp(X)
densityVoronoi(Y)
}