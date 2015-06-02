\name{summary.gcmr}

\alias{print.gcmr}
\alias{summary.gcmr}
\alias{print.summary.gcmr}

\alias{coef.gcmr}
\alias{vcov.gcmr}
\alias{bread.gcmr}
\alias{estfun.gcmr}

\alias{coeftest.gcmr}
\alias{logLik.gcmr}
\alias{terms.gcmr}
\alias{model.frame.gcmr}
\alias{model.matrix.gcmr}

\title{Methods for gcmr Objects}

\description{
  Methods for extracting information from fitted beta
  regression model objects of class \code{"gcmr"}.
}

\usage{
\method{summary}{gcmr}(object, \dots)

\method{coef}{gcmr}(object, \dots)
\method{vcov}{gcmr}(object, \dots)
\method{bread}{gcmr}(x, \dots)
\method{estfun}{gcmr}(x, \dots)
}
\arguments{
\item{object, x}{a fitted marginal regression model of class \code{\link{gcmr}}.}
\item{\dots}{additional arguments, but currently not used.}
}
\value{
The function \code{summary.gcmr} returns an object of class "summary.glm", a list with some components of the \code{gcmr} object, plus 
\item{coefficients}{a list with components \code{marginal} and \code{copula} containing the maximum likelihood estimates of the marginal and Gaussian copula parameters, respectively.}
\item{aic}{Akaike Information Criterion.}

Function \code{\link{coef}} returns the estimated coefficients and \code{\link{vcov}} their variance-covariance matrix. Functions \code{\link{bread}} and \code{\link{estfun}} extract the components of the robust sandwich variance matrix that can be computed with the \code{\link{sandwich}} package (Zeileis, 2004; 2006).
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
\examples{
data(epilepsy)
fit <- gcmr(counts ~ offset(log(time)) + visit + trt + visit:trt, data = epilepsy, 
subset = (id != 49), marginal = negbin.marg, cormat = cluster.cormat(id, "ar1"), 
options=gcmr.options(seed=123, nrep=c(25,100) ))
summary(fit)
}