\name{Eastbourne}
\alias{Eastbourne}
\docType{data}
\title{
Eastbourne traffic accident data
}
\description{
This dataset represents the spatio-temporal locations of traffic accidents in
the down-town of Eastbourne (UK) in the period of 2005-2010. The network was provided by ``OS OpenData" at www.ordnancesurvey.co.uk and is usable under the terms of the OS OpenData license. The traffic locations were collected by the UK Department for Transport at www.data.gov.uk and obtained through kaggle at www.kaggle.com.

The dataset \code{\link{Eastbourne}} is an object of class \code{\link{stlpp}}.

}
\usage{
data(Eastbourne)
}
\source{
Usability: The network of Eastbourne was provided by OS OpenData and contains OS data © Crown copyright and database right (2018). The traffic accident locations in Eastbourne were collected by the UK Department for Transport and were provided by kaggle.

This data is a part of enitre data which is selected and converted to this format by Mehdi Moradi.
}
\references{
Moradi, M.M. and Mateu, J. (2019). First and second-order characteristics of spatio-temporal point processes on linear networks. Journal of Computational and Graphical Statistics. In press.
}

\author{
Mehdi Moradi <m2.moradi@yahoo.com>
}
\seealso{
\code{\link{stlpp}}
}
\examples{
data(Eastbourne)
plot(Eastbourne)
}