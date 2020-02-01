\name{STLg}
\alias{STLg}
\title{Pair correlation function for spatio-temporal point processes on linear networks}
\usage{
STLg(X,r=NULL,t=NULL,nxy=10)
}
\description{Pair correlation function for spatio-temporal point processes on linear networks.}
\arguments{
  \item{X}{a realisation of a spatio-temporal point processes on a linear networks. }
  \item{r}{values of argument r where pair correlation function will be evaluated. optional.}
  \item{t}{values of argument t where pair correlation function will be evaluated. optional.}
    \item{nxy}{pixel array dimensions. optional.}

}

\details{
This function calculates the pair correlation function for a homogeneous spatio-temporal point processes on a linear network.
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
\code{\link{pcf}}, \code{\link{STLK}}
}
\examples{
X <- rpoistlpp(.2,a=0,b=5,L=easynet)
g <- STLg(X)
plot(g)
}
