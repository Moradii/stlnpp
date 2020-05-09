\name{as.stlpp}
\alias{as.stlpp}
\title{Convert data to a space-time point pattern on a linear network}
\usage{
as.stlpp(x,y,t,L)
}
\description{Convert data to a space-time point pattern on a linear network}

\arguments{
  \item{x,y,t}{vectors of cartesian coordinates and time occurance. Alternatively, x can be of classes \code{\link{data.frame}}, \code{\link{ppp}} and \code{\link{lpp}}.}
  \item{L}{linear network (object of class \code{\link{linnet}}) on which the points lie.}
}
\details{
This function converts data to an object of class stlpp.

Data can be of formats:

\itemize{
\item x is of class class \code{\link{data.frame}} with three columns. Then columns are considered as cartesian coordinates (i.e. x,y,t) and they will be converted to a spatio-temporal point pattern on the linear network L.

\item x is a planar point pattern (class \code{\link{ppp}}). Then x will be converted to a spatio-temporal point pattern on the linear network L and with coresponding time vector t.

\item x is a linear point pattern (class \code{\link{lpp}}). Then x will be converted to a spatio-temporal point pattern on the linear network L and with coresponding time vector t.

\item x,y,t are vectors of same length where x,y are living on the corresponding network L.
}
}

\value{
A spatio-temporal point pattern on a linear network. An object of class  \code{\link{stlpp}}.
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}

\seealso{
\code{\link{as.lpp.stlpp}}, \code{\link{runifpointOnLines}}, \code{\link{as.lpp}}
}
\examples{
data(easynet)
x <- runifpointOnLines(40, easynet)
t1 <- sample(1:10,40,replace=TRUE)
Y <- as.stlpp(x,t=t1,L=easynet)

Z <- as.lpp.stlpp(Y)
t2 <- sample(1:10,40,replace=TRUE)
W <- as.stlpp(Z,t=t2)
}
