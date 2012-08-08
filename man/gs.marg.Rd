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
\code{gaussian.marg} \tab Gaussian margins. \cr 
\code{binomial.marg} \tab binomial margins. \cr
\code{poisson.marg} \tab Poisson margins. \cr 
\code{negbin.marg} \tab negative binomial margins. \cr
\code{Gamma.marg} \tab Gamma margins. \cr
\code{weibull.marg} \tab weibull margins. \cr
}
}
\usage{
gaussian.marg(link = "identity")
binomial.marg(link = "logit")
poisson.marg(link = "log")
negbin.marg(link = "log")
Gamma.marg(link = "inverse")
weibull.marg(link = "log")
}

\arguments{
  \item{link}{
a specification for the model link function. See \code{\link{family}} for the special case of generalized linear models. 
}
}
\details{
For binomial marginals specified by \code{bn.marg} the response is specified as a factor when the first level denotes failure and all others success or as a two-column matrix with the columns giving the numbers of successes and failures.

Negative binomial margins implemented in \code{nb.marg} are parameterized such that \eqn{var(Y)=E(Y)+k E(Y)^2}. 
For back-compatibility with previous versions of this package, short names for the margins \code{gs.marg}, \code{bn.marg}, \code{ps.marg} and \code{nb.marg} remain valid as an alternative to (preferred) longer versions  \code{gaussian.marg}, \code{binomial.marg}, \code{poisson.marg} and \code{negbin.marg}.
}
\value{
An object of class \code{\link{marginal.gcmr}} representing the marginal component. 
}
\references{
Masarotto, G. and Varin, C. (2011). Gaussian copula marginal regression. Preprint. 
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}
}
\keyword{regression}
\keyword{nonlinear}
 