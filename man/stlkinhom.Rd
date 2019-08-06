\name{STLkinhom}
\alias{STLkinhom}
\title{spatio-temporal inhomogeneous K-function for point processes on linear networks}
\usage{
STLkinhom(X,lambda=lambda,normaliz=FALSE,seqr=NULL,seqt=NULL)
}
\description{
spatio-temporal inhomogeneous K-function for point processes on linear networks
}
\arguments{
  \item{X}{a realisation of a spatio-temporal point processes on a linear networks. }
  \item{lambda}{values of estimated intensity.}
  \item{normalize}{normalization factor to be considered.}
  \item{seqr}{values of argument r where pair correlation function will be evaluated. optional.}
  \item{seqt}{values of argument t where pair correlation function will be evaluated. optional.}
}
\details{
This function calculates the inhomogeneous K-function for a spatio-temporal point processes on a linear network.
}

\value{
an object of class \code{\link{sumstlpp}}.
}
\references{
Moradi, M.M. and Mateu, J. (2019). First and second-order characteristics of spatio-temporal point processes on linear networks.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{STLg}}, \code{\link{STLK}}, \code{\link{STLginhom}}
}
\examples{
X = rpoistlpp(.2,a=0,b=5,L=easynet)
lambda = density(X)
STLKinhom(X,lambda=lambda)
}

