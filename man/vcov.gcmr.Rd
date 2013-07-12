\name{vcov.gcmr}
\alias{vcov.gcmr}
\title{
Variance-Covariance Matrix of Parameters for a Fitted Gaussian Copula Marginal Regression Model
}
\description{
This method extracts or computes the variance-covariance matrix for parameters in a fitted marginal regression model.
}
\usage{
\method{vcov}{gcmr}(object, type = c("hessian", "sandwich", "vscore", "cluster", "hac"), ...)
}
\arguments{
\item{object}{
a fitted marginal regression model of class \code{\link{gcmr}}.
}
\item{type}{
type of variance-covariance matrix. This can be one of the following:
\tabular{ll}{
\code{hessian} \tab inverse of the observed Fisher information (default). \cr
\code{sandwich}\tab sandwich variance matrix. \cr
\code{vscore} \tab cross-product of scores. \cr
\code{cluster} \tab sandwich variance matrix for cluster/longitudinal data. \cr
\code{hac} \tab heteroskedasticity and autocorrelation consistent estimator for time series of Andrews (1991).
}
}
  \item{\dots}{
further arguments passed to function \code{\link{vcovHAC}}.
}
}
\details{
\code{sandwich} and \code{vscore} are based on the predictive decomposition of the joint density. \code{cluster} uses the decomposition of the data in independent subject-specific blocks. \code{hac} is appropriate for time series and uses the implementation in package \code{\link{sandwich}} (Zeileis, 2006).
}
\value{
A matrix containing the estimated covariances between the parameter estimates of a fitted \code{\link{gcmr}} model.
}
\references{
Andrews, D.W.K. (1991). Heteroskedasticity and autocorrelation consistent covariance matrix estimation. \emph{Econometrica}, \bold{59}, 817--858.

Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics}, \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.

Zeleis, A. (2006). Object-oriented computation of sandwich estimators. \emph{Journal of Statistical Software} \bold{16}, issue 9.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}, \code{\link{sandwich}}, \code{\link{vcovHAC}}
}
\examples{
data(polio)
## working independence       
m <- gcmr( y ~ ., data = polio, marginal = nb.marg(), cormat = ind.cormat() )
m
## HAC variance-covariance matrix
round( vcov( m, "hac" ), 2 )
}
\keyword{regression}
\keyword{nonliner}
