\name{plot.gcmr}
\alias{plot.gcmr}

\title{Plot Diagnostics for Gaussian Copula Marginal Regression}

\description{
  Various types of diagnostic plots for Gaussian copula regression.
}

\usage{
\method{plot}{gcmr}(x, which = if (!time.series) 1:4 else c(1, 3, 5, 6), 
    caption = c("Residuals vs indices of obs.", "Residuals vs linear predictor", 
        "Normal plot of residuals", "Predicted vs observed values", 
        "Autocorrelation plot of residuals", "Partial ACF plot of residuals"), 
    	main = "", ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
    	level = 0.95, col.lines = "gray",
	time.series = inherits(x$cormat, "arma.gcmr"), \ldots)
}
\arguments{
\item{x}{a fitted model object of class \code{\link{gcmr}}.}
\item{which}{select one, or more, of the six available plots. The default choice adapts to the correlation structure and selects four plots depending on the fact that the data are a regular time series or not.}
\item{caption}{captions to appear above the plots.}
\item{main}{title to each plot in addition to the above caption.}
\item{ask}{if \code{TRUE}, then the user is asked before each plot.}
\item{level}{confidence level in the normal probability plot. The default is \code{0.95}.}
\item{col.lines}{color for lines. The default is \code{"gray"}.}
\item{time.series}{if \code{TRUE}, four plots suitable for time series data are displayed. The default is \code{TRUE} when the correlation matrix corresponds to that of ARMA(p,q) process and \code{FALSE} otherwise.}
\item{\dots}{other parameters to be passed through to plotting functions.}
}
\details{
The plot method for \code{\link{gcmr}} objects produces six types of diagnostic plots selectable through the \code{which} argument. Available choices are: Quantile residuals vs indices of the observations (\code{which=1}); Quantile residuals vs linear predictor (\code{which=2}); Normal probability plot of quantile residuals (\code{which=3}); Fitted vs observed values (\code{which=4}); Autocorrelation plot of quantile residuals (\code{which=5}); Partial autocorrelation plot of quantile residuals (\code{which=6}). The latter two plots make sense for regular time series data only.

The normal probability plot is computed via function \code{\link{qqPlot}} from the package \code{\link{car}} (Fox and Weisberg, 2011).
}
\references{
Fox, J. and Weisberg, S. (2011). \emph{An {R} Companion to Applied Regression}. Second Edition. Thousand Oaks CA: Sage. \url{http://socserv.socsci.mcmaster.ca/jfox/Books/Companion}.

Masarotto, G. and Varin, C. (2012). Gaussian copula marginal regression. \emph{Electronic Journal of Statistics} \bold{6}, 1517--1549. \url{http://projecteuclid.org/euclid.ejs/1346421603}.

Masarotto, G. and Varin C. (2017). Gaussian Copula Regression in R. \emph{Journal of Statistical Software}, \bold{77}(8), 1--26. \doi{10.18637/jss.v077.i08}.
}
\author{
Guido Masarotto and Cristiano Varin.
}
\seealso{
\code{\link{gcmr}}.
}
\examples{
## beta regression with ARMA(1,3) errors
data(HUR)
trend <- scale(time(HUR))
m <- gcmr(HUR ~ trend | trend, marginal = beta.marg, cormat = arma.cormat(1, 3))
## normal probability plot
plot(m,  3)
## autocorrelation function of residuals
plot(m,  5)
}