\name{se}
\alias{se}
\title{
Standard Errors of a Fitted Gaussian Copula Marginal Regression Model
}
\description{
Returns the standard errors vector of a fitted Gaussian copula marginal regression model. For a list of possible options see \code{\link{vcov.gcmr}}. 
}
\usage{
se(x, ...)
}
\arguments{
  \item{x}{
a fitted Gaussian copula marginal regression model of class \code{\link{gcmr}}.
}
  \item{\dots}{
further arguments passed to function \code{\link{vcov.gcmr}}.
}
}
\value{
A vector containing the estimated standard errors of the parameters of a fitted \code{\link{gcmr}} model.
}
\references{
Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}, \code{\link{vcov.gcmr}}.
}
\examples{
## beta regression with ARMA(1,3) errors
data(HUR)
trend <- scale(time(HUR))
m <- gcmr(HUR ~ trend | trend, marginal = beta.marg, 
cormat = arma.cormat(1, 3))
## various type of standard errors
res <- cbind( se( m, "he" ), se( m, "v" ), se( m, "s" ), se( m, "hac" ) )
colnames(res) <- c( "hessian", "vscore", "sandwich", "HAC" )
res
}
\keyword{regression}
\keyword{nonliner}
