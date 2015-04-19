\name{estfun.gcmr}
\alias{estfun.gcmr}
\alias{bread.gcmr}
\title{
Components for Sandwich Variance Estimation in Gaussian Marginal Regression Models
}
\description{
Extracts components for computation of the sandwich variance-covariance matrix  (Zeileis, 2004; 2006). 
}
\usage{
\method{estfun}{gcmr}(x, ...)
\method{bread}{gcmr}(x, ...)
}
\arguments{

  \item{x}{
a fitted marginal regression model of class \code{\link{gcmr}}.
}
  \item{\dots}{
further arguments passed to methods.
}
}
\details{
These functions are passed to \code{\link{sandwich}} (Zeileis, 2004; 2006).
}
\value{
\code{estfun}: a matrix containing the empirical estimating functions. 

\code{bread}: a matrix containing the `bread' of the sandwich estimator. 
}
\references{
Zeileis, A. (2004). Econometric computing with HC and HAC covariance matrix estimators. \emph{Journal of Statistical Software} \bold{11}, issue 10.

Zeileis, A. (2006). Object-oriented computation of sandwich estimators. \emph{Journal of Statistical Software} \bold{16}, issue 9.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{bread}}, \code{\link{estfun}}, \code{\link{gcmr}}, \code{\link{sandwich}}.
}
\keyword{regression}
\keyword{nonlinear}
