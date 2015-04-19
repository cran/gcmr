\name{residuals.gcmr}
\alias{residuals.gcmr}
\title{Quantile Residuals for Gaussian Copula Marginal Regression}
\usage{
\method{residuals}{gcmr}(object, type=c("conditional","marginal"),
          method=c("random","mid"),...)
}
\arguments{
  \item{object}{an object of class \code{gcmr}, typically the result of
    a call to \code{\link{gcmr}}.}
  \item{type }{the type of quantile residuals which should be returned.
    The alternatives are: \code{"conditional"} (default) and \code{"marginal"}.}
  \item{method }{different methods available for quantile residuals in case of discrete responses:
    \code{"random"} for randomized quantile residuals (default), and 
    \code{"mid"} for mid interval quantile residuals as defined in Zucchini and MacDonald (2009).}
\item{...}{further arguments passed to or from other methods.}
}
\description{
  Computes various type of quantile residuals for validation of a fitted Gaussian copula marginal regression model described in Masarotto and Varin (2012).
}
\details{
 Quantile residuals are defined in Dunn and Smyth (1996). Two different types are available:
\tabular{ll}{
\code{conditional} \tab quantile residuals that account for the dependence. \cr
\code{marginal} \tab quantile residuals that do not account for the dependence.
}
Conditional quantile residuals are normal quantiles of Rosenblatt (1952) transformations and they are appropriate for validation of the marginal regression models discussed by Masarotto and Varin (2012). If the responses are discrete, then the conditional quantile residuals are not well defined. This difficulty is overcame by randomized quantile residuals available through option \code{method="random"}. Alternatively, Zucchini and MacDonald (2009) suggest the use of mid interval quantile residuals (\code{method="mid"}). 
}
\seealso{
  \code{\link{gcmr}}
}
\references{
Dunn, P.K. and Smyth, G.K. (1996). Randomized quantile residuals. \emph{Journal of Computational and Graphical Statistics} \bold{5}, 236--244.

Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.

Rosenblatt, M. (1952). Remarks on a multivariate transformation. \emph{The Annals of Mathematical Statistics} \bold{23}, 470--472.

Zucchini, W. and MacDonald, I.L. (2009). \emph{Hidden Markov Models for Time Series}. Chapman and Hall/CRC.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\note{
Differently from randomized quantile residuals, mid quantile residuals are \bold{not} realizations of incorrelated standard normal variables under model conditions.

It is appropriate to inspect several sets of randomized quantile residuals before to take a decision about the model.

See Masarotto and Varin (2012) for more details.
}
\examples{
## spatial binomial data
data(malaria)
D <- sp::spDists(cbind(malaria$x, malaria$y))/1000
m <- gcmr(cbind(cases, size-cases) ~ netuse+I(green/100)+phc, data=malaria, 
marginal=binomial.marg, cormat=matern.cormat(D))
res <- residuals(m)
## normal probability plot
qqnorm(res)
qqline(res)
## or better via plot.gcmr
plot(m, which = 3)
}
\keyword{models}
\keyword{regression}

