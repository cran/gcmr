\name{cormat.gcmr}
\title{Correlation Matrices for Gaussian Copula Regression Models}
\alias{cormat.gcmr}
\description{
  Class of correlation matrices available in the \code{\link{gcmr}} package. 
}
\value{
At the moment, the following are implemented:
\tabular{ll}{ 
\code{\link{ind.cormat}} \tab working independence. \cr
\code{\link{arma.cormat}} \tab ARMA(p,q). \cr 
\code{\link{cluster.cormat}} \tab longitudinal/clustered data. \cr
\code{\link{matern.cormat}} \tab  Matern spatial correlation. \cr
}
}
\references{
Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{\code{\link{gcmr}},
  \code{\link{ind.cormat}},
  \code{\link{arma.cormat}},
  \code{\link{cluster.cormat}}, 
  \code{\link{matern.cormat}}.
}
\keyword{regression}
\keyword{nonlinear}
