\name{rpoistpp}
\alias{rpoistpp}
\title{
Simulating temporal poisson point patterns
}
\description{
simulating realisations of a temporal poisson point process.}
\usage{
rpoistpp(lambda,a,b,check=FALSE,lmax=NULL,nsim=1)
}
\arguments{
\item{lambda}{
intensity of the point process. It can be either a number, a function of location and time, or an object of class \code{tppint}.
}
\item{a}{
lower bound of time period.
}
\item{b}{
upper bound of time period.
}
\item{check}{
Logical value indicating whether to check that all the points lie inside the specified time period.
}
\item{lmax}{
upper bound for the values of labmda. This is optinal.
}
\item{nsim}{
number of simulated patterns to generate.
}
}
\details{
This function generates realisations of a temporal poisson point process based on a given intensity function lambda and lower/upper bounds a and b.
}
\value{
an object of class \code{\link{tpp}} if nsim=1, otherwise a list of objects of class \code{\link{tpp}}.}
\references{
Moradi, M.M. and Mateu, J. (2019). First and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics. In press.}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{rpoistlpp}}
}
\examples{
f <- function(t){0.1*exp(t)}
X <- rpoistpp(f,a=1,b=10)
}