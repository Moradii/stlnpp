\name{rpoistlpp}
\alias{rpoistlpp}
\title{Simulating spatio-temporal poisson point processes on a linear network}
\usage{
rpoistlpp(lambda,a,b,L,check=FALSE,lmax=NULL,nsim=1)
}
\description{
simulating a realisation of a spatio-temporal poisson point process on a linear network.
}

\arguments{
  \item{lambda}{intensity of the point process. It can be either a number, function of location and time, or an abject of class \code{stlppint}.}
  \item{a}{lower bound of time period.}
  \item{b}{upper bound of time period.}
  \item{L}{a linear network.}
  \item{check}{Logical value indicating whether to check that all the (x,y) points lie inside the specified window. see \code{\link{ppp}}.}
  \item{lmax}{upper bound for the values of labmda. This is optinal.}
  \item{nsim}{number of simulated patterns to generate.}
}

\details{
This function generates realisations of a spatio-temporal poisson point process on a linear network based on an intensity function lambda and lower/upper bounds a and b.
}

\value{
an object of class \code{\link{stlpp}} if nsim=1, otherwise a list of objects of class \code{\link{stlpp}}.
}
\references{
Moradi, M.M. and Mateu, J. (2019). First and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics. In press.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{density.stlpp}}
}
\examples{
X <- rpoistlpp(.2,a=0,b=5,L=easynet)
X
}

