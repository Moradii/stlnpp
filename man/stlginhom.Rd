\name{STLginhom}
\alias{STLginhom}
\title{spatio-temporal inhomogeneous pair correlation function for point processes on linear networks}
\usage{
STLginhom(X,lambda,normalize=FALSE,seqr=NULL,seqt=NULL)
}
\description{
spatio-temporal inhomogeneous pair correlation function for point processes on linear networks.
}
\arguments{
  \item{X}{a realisation of a spatio-temporal point processes on a linear networks. }
  \item{lambda}{values of estimated intensity.}
  \item{normalize}{normalization factor to be considered.}
  \item{seqr}{values of argument r where pair correlation function will be evaluated. optional.}
  \item{seqt}{values of argument t where pair correlation function will be evaluated. optional.}
}

\details{
This function calculates the inhomogeneous pair correlation function for a spatio-temporal point processes on a linear network.
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
\code{\link{STLg}}, \code{\link{STLK}}, \code{\link{STLKinhom}}
}
\examples{
X = rpoistlpp(.2,a=0,b=5,L=easynet)
d = density(X)
STLginhom(X,lambda=d)
}
