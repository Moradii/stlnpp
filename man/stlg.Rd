\name{STLg}
\alias{STLg}
\title{spatio-temporal pair correlation function for point processes on linear networks}
\usage{
STLg(X,seqr=NULL,seqt=NULL)
}
\description{spatio-temporal pair correlation function for point processes on linear networks.}
\arguments{
  \item{X}{a realisation of a spatio-temporal point processes on a linear networks. }
  \item{seqr}{values of argument r where pair correlation function will be evaluated. optional.}
  \item{seqt}{values of argument t where pair correlation function will be evaluated. optional.}
}

\details{
This function calculates the pair correlation function for a homogeneous spatio-temporal point processes on a linear network.
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
\code{\link{pcf}}, \code{\link{STLK}}
}
\examples{
X = rpoistlpp(.2,a=0,b=5,L=easynet)
STLg(X)
}
