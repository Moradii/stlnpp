\name{density.stlpp}
\alias{density.stlpp}
\title{Kernel estimation of intensity of space-time point patterns on a linear network}
\usage{
\method{density}{stlpp}(X,lbw,tbw,at=c("points","pixels"),dimt=512,...)
}
\description{
Kernel density estimation of a spatio-temporal point pattern on a linear network.
}

\arguments{
  \item{X}{an object of class \code{\link{stlpp}}}
  \item{lbw}{netwrok smoothing bandwidth}
  \item{tbw}{time smoothing bandwidth}
  \item{at}{string specifying whether to compute the intensity values at a grid of pixel locations and time (at="pixels") or only at the points of x (at="points"). default is to estimate the intensity at pixels}
  \item{dimt}{the number of equally spaced points at which the temporal density is to be estimated. see \link[stats]{density}}
  \item{...}{arguments passed to \code{\link{density.lpp}}}
}

\details{
Kernel smoothing is applied to the spatio-temporal point pattern X using methods in Moradi et al (2019). If lbw and tbw are not given, then they will be selected using \code{\link{bw.nrd0}} and \code{\link{bw.scott.iso}} respectively.
}

\value{
if \code{at="points"}: a vector of intensity values at the data points of X.

if \code{at="pixels"}: a list of images on linear network. Each image represents an estimated saptio-temporal intensity at a fixed time.
check the attributes for more accommodated outputs.
}

\references{
Moradi, M.M. and Mateu, J. (2019). First and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics. In press.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{density}}, \code{\link{density.lpp}}, \code{\link{bw.nrd0}}, \code{\link{bw.scott.iso}}
}
\examples{
X <- rpoistlpp(.2,a=0,b=5,L=easynet)
density(X)
}
