\name{gaussian.marg}
\alias{gaussian.marg}
\alias{binomial.marg}
\alias{poisson.marg}
\alias{negbin.marg}
\alias{Gamma.marg}
\alias{weibull.marg}
\alias{gs.marg}
\alias{bn.marg}
\alias{ps.marg}
\alias{nb.marg}
\title{
Sets the Marginals in Gaussian Copula Marginal Regression Models 
}
\description{
These functions set the marginals in Gaussian copula marginal regression models.

At the moment, the following are implemented:\tabular{ll}{ 
\code{binomial.marg} \tab binomial margins. \cr
\code{Gamma.marg} \tab Gamma margins. \cr
\code{gaussian.marg} \tab Gaussian margins. \cr 
\code{negbin.marg} \tab negative binomial margins. \cr
\code{poisson.marg} \tab Poisson margins. \cr 
\code{weibull.marg} \tab Weibull margins. \cr
}
}
\usage{
binomial.marg(link = "logit")
Gamma.marg(link = "inverse")
gaussian.marg(link = "identity")
negbin.marg(link = "log")
poisson.marg(link = "log")
weibull.marg(link = "log")
}

\arguments{
  \item{link}{
a specification for the model link function. See \code{\link{family}} for the special case of generalized linear models. 
}
}
\details{
For binomial marginals specified by \code{binomial.marg} the response is specified as a factor when the first level denotes failure and all others success or as a two-column matrix with the columns giving the numbers of successes and failures.

Negative binomial margins implemented in \code{negbin.marg} are parameterized such that \eqn{var(Y)=E(Y)+k E(Y)^2}. 

For back-compatibility with previous versions of this package, short names for the margins \code{bn.marg}, \code{gs.marg}, \code{nb.marg}, and \code{ps.marg} remain valid as an alternative to (preferred) longer versions \code{binomial.marg},  \code{gaussian.marg}, \code{negbin.marg}, and \code{poisson.marg}.
}
\value{
An object of class \code{\link{marginal.gcmr}} representing the marginal component. 
}
\references{
Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics}, \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}
}
\keyword{regression}
\keyword{nonlinear}
 