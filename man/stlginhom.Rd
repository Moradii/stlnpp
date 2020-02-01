\name{STLginhom}
\alias{STLginhom}
\title{Inhomogeneous pair correlation function for spatio-temporal point processes on linear networks}
\usage{
STLginhom(X,lambda,normalize=FALSE,r=NULL,t=NULL,nxy=10)
}
\description{
Inhomogeneous pair correlation function for spatio-temporal  point processes on linear networks.
}
\arguments{
  \item{X}{a realisation of a spatio-temporal point processes on a linear networks. }
  \item{lambda}{values of estimated intensity at data points.}
  \item{normalize}{normalization factor to be considered.}
  \item{r}{values of argument r where pair correlation function will be evaluated. optional.}
  \item{t}{values of argument t where pair correlation function will be evaluated. optional.}
    \item{nxy}{pixel array dimensions. optional.}

}

\details{
This function calculates the inhomogeneous pair correlation function for a spatio-temporal point processes on a linear network.
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
\code{\link{STLg}}, \code{\link{STLK}}, \code{\link{STLKinhom}}
}
\examples{
X <- rpoistlpp(.2,a=0,b=5,L=easynet)
d <- density(X,at="points")
g <- STLginhom(X,lambda=d,normalize=TRUE)
plot(g)
}
