\name{tpp}
\alias{tpp}
\title{
Create a temporal point pattern
}
\description{
Create an object of class \code{\link{tpp}} that represents a temporal point pattern.}
\usage{
tpp(X,a,b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{
an object of class \code{\link{numeric}}, \code{\link{integer}} or \code{\link{vector}}
}
\item{a}{
lower band
}
\item{b}{
upper bound
}
}
\details{
Create a temporal point pattern.}
\value{
An object of class tpp.
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{stlpp}}
}

\examples{
tpp(runif(10))
}