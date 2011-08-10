\name{gs.marg}
\alias{gs.marg}
\alias{bn.marg}
\alias{ps.marg}
\alias{nb.marg}
\alias{sn.marg}
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
\code{sn.marg} \tab skew-normal margins (Azzalini, 1985).
}
}
\usage{
gs.marg(link = "identity")
bn.marg(link = "logit")
ps.marg(link = "log")
nb.marg(link = "log")
sn.marg(link = "identity")
}

\arguments{
  \item{link}{
a specification for the model link function. See \code{\link{family}} for the special case of generalized linear models. 
}
}
\details{
Binomial margins implemented in \code{bn.marg} require that y is the form of a two-column matrix with the columns giving the numbers of successes and failures.

Negative binomial margins implemented in \code{nb.marg} are parameterized such that \eqn{var(Y)=E(Y)+k E(Y)^2}. 

Skew-normal margins implemented in \code{sn.marg} are inherited from package \code{\link{SN}} and use the centered parameterization (Azzalini and Capitanio, 1999).

\bold{Warning}. In order to stabilize numerical optimization, constrained parameters are transformed. The scale parameter sigma of normal and skew-normal margins is log-transformed, while the shape parameter of skew-normal margins is transformed in the logit scale \eqn{log((0.995-shape)/(shape-0.995))}.
}
\value{
An object of class \code{\link{marginal.gcmr}} representing the marginal component. 
}
\references{
Azzalini, A. (1985). A class of distributions which includes the normal ones. \emph{Scandinavian Journal of Statistics} \bold{12}, 171--178.

Azzalini, A. and Capitanio, A. (1999). Statistical applications of the multivariate skew-normal distribution. \emph{Journal of the Royal Statistical Society B} \bold{61}, 579--602.
 
Masarotto, G. and Varin, C. (2011). Gaussian copula marginal regression. Preprint. 
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}, \code{\link{dsn}}
}
\keyword{regression}
\keyword{nonlinear}
 