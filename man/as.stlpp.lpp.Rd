\name{as.lpp.stlpp}
\alias{as.lpp.stlpp}

\title{
Methods for space-time point patterns on a linear network.
}
\description{
This function projects an object of class \code{\link{stlpp}} to a linear network.
}

\usage{
\method{as.lpp}{stlpp}(x,...)
}

\arguments{
\item{x}{an object of class \code{\link{stlpp}}
}
\item{...}{
arguments passed to \code{\link{as.lpp}}
}
}
\details{
This function projects the space-time point pattern x on linear network L into L, giving its corresponding spatial point pattern.
}
\value{
An object of class \code{\link{lpp}}.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}



\seealso{
\code{\link{as.stlpp}}, \code{\link{lpp}}, \code{\link{as.lpp}}
}

\examples{
data(easynet)
x <- runifpointOnLines(40, easynet)
t1 <- sample(1:10,40,replace=TRUE)
Y <- as.stlpp(x,t=t1,L=easynet)
as.lpp.stlpp(Y)
}
