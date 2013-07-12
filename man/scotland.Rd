\name{scotland}
\alias{scotland}
\docType{data}
\title{
Scotland Lip Cancer Data
}
\description{
Male lip cancer in Scotland counties between 1975-1980.
}
\usage{data(scotland)}
\format{
A data frame with the 56 observations with the following variables
\tabular{ll}{

\code{observed} \tab observed cases in each county. \cr
\code{expected} \tab excpected cases in each county. \cr
\code{AFF} \tab proportion of the population employed in agriculture, fishing, or forestry. \cr
\code{latitude} \tab county latitude. \cr
\code{longitude} \tab county longitude. \cr
}
}
\source{
Waller, L.A. and Gotway, C.A. (2004) \emph{Applied Spatial Statistics for Public Health Data}. New York: John Wiley and Sons.
\url{http://www.sph.emory.edu/~lwaller/book/ch2/scotland.dat}
}
\references{
Clayton D. and Kaldor J. (1987) Empirical Bayes estimates of age-standardized relative risks for use in disease mapping. \emph{Biometrics} \bold{43}, 671--681.
}
\examples{
data(scotland)
}
\keyword{datasets}
