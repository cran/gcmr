\name{gaussian.marg}
\alias{gaussian.marg}
\alias{beta.marg}
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
Marginals in Gaussian Copula Marginal Regression Models 
}
\description{
These functions set the marginals in Gaussian copula marginal regression models. 
}
\usage{
beta.marg(link = "logit")
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
Beta marginals specified by \code{beta.marg} are parametrized in terms of mean and dispersion as in \code{\link{betareg}}. See Cribari-Neto and Zeileis (2010) and Ferrari and Cribari-Neto (2004). 

For binomial marginals specified by \code{binomial.marg}, the response is specified as a factor when the first level denotes failure and all others success or as a two-column matrix with the columns giving the numbers of successes and failures.

Negative binomial marginals implemented in \code{negbin.marg} are parametrized such that \eqn{var(Y)=E(Y)+k E(Y)^2}. 

For back-compatibility with previous versions of the \code{gcmr} package, short names for the marginals \code{bn.marg}, \code{gs.marg}, \code{nb.marg}, and \code{ps.marg} remain valid as an alternative to (preferred) longer versions \code{binomial.marg},  \code{gaussian.marg}, \code{negbin.marg}, and \code{poisson.marg}.
}
\value{
An object of class \code{\link{marginal.gcmr}} representing the marginal component. 
}
\references{
Cribari-Neto, F. and Zeileis, A. (2010). Beta regression in R. \emph{Journal of Statistical Software} \bold{34}, 1--24. \url{http://www.jstatsoft.org/v34/i02/}.

Ferrari, S.L.P. and Cribari-Neto, F. (2004). Beta regression for modeling rates and proportions. \emph{Journal of Applied Statistics} \bold{31} (7), 799--815.

Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.

Masarotto, G. and Varin C. (2017). Gaussian Copula Regression in R. \emph{Journal of Statistical Software}, \bold{77}(8), 1--26. \doi{10.18637/jss.v077.i08}.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}, \code{\link{betareg}}.
}
\keyword{regression}
\keyword{nonlinear}
 