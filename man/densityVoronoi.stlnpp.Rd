\name{densityVoronoi.stlpp}
\alias{densityVoronoi.stlpp}
\title{
Intensity estimate of spatio-temporal point pattern using Voronoi-Dirichlet tessellation}
\description{
adaptive intensity estimation for spatio-temporal point patterns on linear networks using Voronoi-Dirichlet tessellation.}
\usage{
\method{densityVoronoi}{stlpp}(X, f = 1, nrep = 1,
                      separable=FALSE,at=c("points","pixels"),
                      dimt=128,...)
}
\arguments{
\item{X}{
an object of class \code{\link{stlpp}}
}
\item{f}{
fraction (between 0 and 1 inclusive) of the data points that will be used to build a tessellation for the intensity estimate
}
\item{nrep}{
number of independent repetitions of the randomised procedure
}
\item{separable}{
logical. If FALSE, it then calculates a pseudo-separable estimate
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
This function computes intensity estimates for spatio-temporal point patterns on linear networks using Voronoi-Dirichlet tessellation. Both first-order separability and pseudo-separability assumptions are accommodated in the function.

If separable=TRUE, the estimated intensities will be a product of the estimated intensities on the network and those on time. Each will be obtained using \code{\link{densityVoronoi.lpp}}. If f=1, the function calculates the estimations based on the original Voronoi intensity estimator.


If separable=FALSE, the estimated intensities will be calculated based on a subsampling technique explained in Mateu et al. (2019). nrep subsamples will be obtained from X based on a given retention probability f, the function \code{\link{densityVoronoi.stlpp}}, considering separable=TRUE and f=1, will be applied to each obtained subsample, and finally, the estimated intensities will be the sum of all obtained estimated intensities from all subsamples divided by the (f * nrep).

}
\value{
If \code{at="points"}: a vector of intensity values at the data points of X.

If \code{at="pixels"}: a list of images on a linear network. Each image represents an estimated saptio-temporal intensity at a fixed time.
}
\references{
Mateu, J., Moradi, M., & Cronie, O. (2019). Spatio-temporal point patterns on linear networks: Pseudo-separable intensity estimation. Spatial Statistics, 100400.
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>,  and Ottmar Cronie.
}

\seealso{
\code{\link{densityVoronoi.lpp}}, \code{\link{density.stlpp}}
}
\examples{
X <- rpoistlpp(0.2,a=0,b=5,L=easynet)
densityVoronoi(X)
}
