## vcov must have class vcov.mr and the following elements:
## - npar: number of parameters. Wow. Really?
## - start(marginal, betastart): computes initial estimates
## - chol(tau,not.na): returns the cholesky factor of the vcov matrix (only for the not.na obs)
## returns NULL if gamma is outside parameter space

## working independence correlation
ind.cormat <- function() {
    ans <- list()
    ans$npar <- 0
    ans$start <- function() double(0)
    ans$chol <- function(tau , not.na) diag(rep(1,sum(not.na)))
    ans$independent <- TRUE
    class( ans ) <- "cormat.gcmr"
    ans
}

## arma(p,q) correlation for time-series
arma.cormat <- function( p , q ) {
    if(p==0 && q==0)
      return( ind.cormat() )
    start <- rep( 0 , p+q )
    names(start) <- c(if ( p ) paste("ar",1:p,sep="") else NULL ,
                      if ( q ) paste("ma",1:q,sep="") else NULL )
    iar <- if ( p ) 1:p else NULL
    ima <- if ( q ) (p+1):(p+q) else NULL
    ans <- list()
    ans$npar <- length(start)
    ans$start <- function() start
    ans$chol <- function( tau , not.na ) {
        if ( ( p && any(Mod(polyroot(c(1,-tau[iar])))<1.01) ) ||
             ( q && any(Mod(polyroot(c(1, tau[ima])))<1.01) ) )
            return( NULL )
        n <- length(not.na)
        rho <- ARMAacf(tau[iar],tau[ima],n-1)
        r <- seq(1,n)[not.na]
        chol(outer( r , r , function(i,j) rho[1+abs(i-j)] ))
    }
    class( ans ) <- "cormat.gcmr"
    ans
}

## clustered data
## assume that it is not possible that all the observations inside a cluster
## can be missing
cluster.cormat <- function(id, type=c("ar1", "ma1", "exch", "unstr")) { 
    type <- match.arg(type)
    if(!length(rle(id)$values)==length(unique(id)))
      stop("length of 'id' must be the same of the number of observations and data must be
sorted in way that observations from the same cluster are contiguous")
    ng <- 1:length(unique(id))
    if (!(length(ng)>1)) stop("only one strata")
    ans <- list(type=type,id=id)
    if(type=="unstr")
        ans$npar <- choose(max(table(id)), 2)
    else
        ans$npar <- 1
    start <- rep(0, ans$npar)
    names(start) <- paste("tau", 1:ans$npar, sep="")
    data <- data.frame(id=id)
    fn <- switch(type, "ar1"=function(g) nlme::corAR1(g, form= ~1|id),
                 "ma1"=function(g) nlme::corARMA(g, form= ~1|id, p=0, q=1),
                 "exch"=function(g) nlme::corCompSymm(g, form= ~1|id),
                 "unstr"=function(g) nlme::corSymm(g, form= ~1|id))
    ans$start <- function() start
    ans$chol <- function(tau, not.na) {
        q <- try(nlme::corMatrix(nlme::Initialize(fn(tau),data=data)),silent=TRUE)
        if (inherits(q,"try-error")) return(NULL)
        g <- split(not.na,id)
        q <- try(lapply(ng,function(i) chol(q[[i]][g[[i]],g[[i]]])),silent=TRUE)
        if (inherits(q,"try-error") ) NULL else q
    }
    class( ans ) <- "cormat.gcmr"
    ans
}

## Matern correlation for spatial data
## D is a distance matrix, alpha is the smoothing parameter
matern.cormat <- function(D, alpha=0.5){
  ans <- list()
  ans$npar <- 1
  start <- median(D) 
  names(start) <- c("tau")
  ans$start <- function( marginal, beta ) start
  ans$chol <- function( tau, not.na ){
    if( any(tau<=0) ) return( NULL )
    S <- geoR:::matern(D, tau, alpha)
    q <- try(chol(S),silent=TRUE)
    if( inherits(q,"try-error") ) NULL else q
  }
  class( ans ) <- "cormat.gcmr"
  ans
}
