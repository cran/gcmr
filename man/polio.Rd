\name{polio}
\alias{polio}
\docType{data}
\title{
Polio Time Series 
}
\description{
Time series of Polio incidences in U.S.A. from 1970 to 1983.
}
\usage{data(polio)}
\format{
A data frame with the 168 monthly observations (from January 1970 to December 1983) with the following variables
\tabular{ll}{
\code{y} \tab time series of polio incidences. \cr
\code{t*10^(-3)} \tab linear trend multiplied by factor \eqn{10^(-3)}. \cr
\code{cos(2*pi*t/12)} \tab cosine annual seasonal component. \cr
\code{sin(2*pi*t/12)} \tab sine annual seasonal component. \cr
\code{cos(2*pi*t/6)} \tab cosine semi-annual seasonal component. \cr
\code{sin(2*pi*t/6)} \tab sine semi-annual seasonal component. \cr
}
}
\source{
Zeger, S.L. (1988). A regression model for time series of counts. \emph{Biometrika} \bold{75}, 822--835.
}
\examples{
data(polio)
}
\keyword{datasets}
