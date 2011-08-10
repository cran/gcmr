\name{marginal.gcmr}
\title{Marginal Component in Marginal Regression Models}
\alias{marginal.gcmr}
\description{
  Class of marginals available in the \code{gcmr} library. 
}
\value{
At the moment, the following are implemented:
\tabular{ll}{ 
\code{\link{gs.marg}} \tab normal margins. \cr 
\code{\link{bn.marg}} \tab binomial margins. \cr
\code{\link{ps.marg}} \tab Poisson margins. \cr 
\code{\link{nb.marg}} \tab negative binomial margins. \cr
\code{\link{sn.marg}} \tab skew-normal margins.  
}
}
\references{
 Masarotto, G. and Varin, C. (2011). Gaussian copula marginal regression. Preprint.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{\code{\link{gcmr}},
  \code{\link{gs.marg}},
  \code{\link{bn.marg}},
  \code{\link{ps.marg}}, 
  \code{\link{nb.marg}},
  \code{\link{sn.marg}}
}
\keyword{regression}
\keyword{nonlinear}
