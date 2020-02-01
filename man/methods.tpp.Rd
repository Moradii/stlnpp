\name{methods.tpp}
\alias{methods.tpp}
\alias{print.tpp}
\alias{print.tppint}
\alias{plot.tpp}
\alias{plot.tppint}
\title{
Methods for temporal point patterns
}
\description{
Methods for temporal point patterns.
}
\usage{
\method{plot}{tpp}(X,xlab="time",ylab="",main = "cumulative number",...)
\method{plot}{tppint}(X,xlab=xlab,xlim=xlim,line=2.5,...)
\method{print}{tpp}(X)
\method{print}{tppint}(X)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
an object of class tpp or tppint.
}
\item{xlab,ylab}{the x,y label of the plot.}
\item{main}{overall title for the plot.}
\item{xlim}{giving the x limits for the plot.}
\item{line}{specifying a value for line overrides the default placement of y label, and places it this many lines outwards from the plot edge.}
\item{...}{graphics parameters passed to plot function.}
}
\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}
\examples{
X <- tpp(sample(c(1:24),200,replace = TRUE))
plot(X)
plot(density(X))
}

