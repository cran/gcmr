\name{marginal.gcmr}
\title{Marginal Component in Marginal Regression Models}
\alias{marginal.gcmr}
\description{
  Class of marginals available in the \code{gcmr} library. 
}
\value{
At the moment, the following are implemented:
\tabular{ll}{ 

\code{\link{beta.marg}} \tab beta margins. \cr
\code{\link{binomial.marg}} \tab binomial margins. \cr
\code{\link{Gamma.marg}} \tab Gamma margins. \cr
\code{\link{gaussian.marg}} \tab normal margins. \cr 
\code{\link{negbin.marg}} \tab negative binomial margins. \cr
\code{\link{poisson.marg}} \tab Poisson margins. \cr 
\code{\link{weibull.marg}} \tab Weibull margins. \cr
}
}
\references{
Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.
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
