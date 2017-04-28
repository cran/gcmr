\name{marginal.gcmr}
\title{Marginals for Gaussian Copula Marginal Regression}
\alias{marginal.gcmr}
\description{
  Class of marginals available in the \code{gcmr} library. 
}
\value{
At the moment, the following are implemented:
\tabular{ll}{ 

\code{\link{beta.marg}} \tab beta marginals. \cr
\code{\link{binomial.marg}} \tab binomial marginals. \cr
\code{\link{Gamma.marg}} \tab Gamma marginals. \cr
\code{\link{gaussian.marg}} \tab Gaussian marginals. \cr 
\code{\link{negbin.marg}} \tab negative binomial marginals. \cr
\code{\link{poisson.marg}} \tab Poisson marginals. \cr 
\code{\link{weibull.marg}} \tab Weibull marginals. \cr
}
}
\references{
Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.

Masarotto, G. and Varin C. (2017). Gaussian Copula Regression in R. \emph{Journal of Statistical Software}, \bold{77}(8), 1--26. \doi{10.18637/jss.v077.i08}.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{\code{\link{gcmr}},
	\code{\link{beta.marg}},		
  	\code{\link{binomial.marg}},
  	\code{\link{gaussian.marg}},
	\code{\link{Gamma.marg}},
  	\code{\link{negbin.marg}},
  	\code{\link{poisson.marg}},   
  	\code{\link{weibull.marg}}.
}
\keyword{regression}
\keyword{nonlinear}
