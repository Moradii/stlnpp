\name{STLK}
\alias{STLK}
\title{K-function for spatio-temporal point processes on linear networks}
\usage{
STLK(X,r=NULL,t=NULL,nxy=10)
}
\description{
K-function for spatio-temporal point processes on linear networks.
}
\arguments{
  \item{X}{a realisation of a spatio-temporal point processes on a linear networks. }
  \item{r}{values of argument r where pair correlation function will be evaluated. optional.}
  \item{t}{values of argument t where pair correlation function will be evaluated. optional.}
    \item{nxy}{pixel array dimensions. optional.}
}

\details{
This function calculates the K-function for a homogeneous spatio-temporal point processes on a linear network.
}

\value{
An object of class \code{sumstlpp}.
}
\references{
Moradi, M.M. and Mateu, J. (2019). First and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics. In press.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{Kest}}, \code{\link{STLg}}
}
\examples{
X <- rpoistlpp(.2,a=0,b=5,L=easynet)
k <- STLK(X)
plot(k)
}

