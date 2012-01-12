\name{gs.marg}
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
\code{gs.marg} \tab normal margins. \cr 
\code{bn.marg} \tab binomial margins. \cr
\code{ps.marg} \tab Poisson margins. \cr 
\code{nb.marg} \tab negative binomial margins. \cr
}
}
\usage{
gs.marg(link = "identity")
bn.marg(link = "logit")
ps.marg(link = "log")
nb.marg(link = "log")
}

\arguments{
  \item{link}{
a specification for the model link function. See \code{\link{family}} for the special case of generalized linear models. 
}
}
\details{
For binomial marginals specified by \code{bn.marg} the response is specified as a factor when the first level denotes failure and all others success or as a two-column matrix with the columns giving the numbers of successes and failures.

Negative binomial margins implemented in \code{nb.marg} are parameterized such that \eqn{var(Y)=E(Y)+k E(Y)^2}. 

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
 