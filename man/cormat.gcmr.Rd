\name{cormat.gcmr}
\title{Correlation Matrices for Gaussian Copula Regression Models}
\alias{cormat.gcmr}
\description{
  Class of correlation matrices available in the \code{gcmr} library. 
}
\value{
At the moment, the following are implemented:
\tabular{ll}{ 
\code{\link{ind.cormat}} \tab working independence. \cr
\code{\link{arma.cormat}} \tab ARMA(p,q).  \cr 
\code{\link{cluster.cormat}} \tab longitudinal/clustered data.  \cr
\code{\link{matern.cormat}} \tab  Matern spatial correlation. \cr
}
}
\references{
 Masarotto, G. and Varin, C. (2011). Gaussian copula marginal regression. Preprint.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{\code{\link{gcmr}},
  \code{\link{ind.cormat}},
  \code{\link{arma.cormat}},
  \code{\link{cluster.cormat}}, 
  \code{\link{matern.cormat}},
}
\keyword{regression}
\keyword{nonlinear}
