\name{methods.tpp}
\alias{methods.tpp}
\alias{print.tpp}
\alias{print.tppint}
\alias{plot.tpp}
\alias{plot.tppint}
\alias{[.tpp}
\alias{[.tppint}
\title{
Methods for temporal point patterns
}
\description{
Methods for temporal point patterns.
}
\usage{
\method{plot}{tpp}(x,xlab="time",ylab="",main = "cumulative number",...)
\method{plot}{tppint}(x,xlab=xlab,xlim=xlim,line=2.5,main="NULL",...)
\method{print}{tpp}(x,...)
\method{print}{tppint}(x,...)
\method{[}{tpp}(x,i)
\method{[}{tppint}(x,i)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
an object of class tpp or tppint.
}
\item{xlab,ylab}{the x,y label of the plot.}
\item{main}{overall title for the plot.}
\item{xlim}{giving the x limits for the plot.}
\item{line}{specifying a value for line overrides the default placement of y label, and places it this many lines outwards from the plot edge.}
\item{i}{
numeric, logical, or an object of class \code{\link{tpp}}
}
\item{...}{graphics parameters passed to plot/print function.}
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}
\examples{
X <- tpp(sample(c(1:24),200,replace = TRUE))
plot(X)
plot(density(X))
}

