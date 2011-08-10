## Marginals

## Marginal must have class marginal.mr and the following elements:
## - start(y,x): compute initial estimates ignoring correlations
## - dp(y,x,beta): evaluate [d,p]
## - q(u,y,x,beta): evaluate quantiles
## - npar(x): number of parameters
## - type: response type (integer or numeric)

# Gaussian
gs.marg <- function(link = "identity" ) {
    fm <- gaussian( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x, offset) {
        m <- glm.fit( x , y, offset=offset, family=fm )
        ans <- c( coef( m ) , log(sqrt( mean( residuals( m )^2 ) ) ))
        names(ans) <- c(dimnames(as.matrix(x))[[2L]],"log.sigma")
        ans
    }
    ans$npar <- function(x) NCOL(x)+1
    ans$dp <- function(y, x, offset, beta) {
        nb <- length(beta)
        mu <- fm$linkinv( x %*% beta[-nb] + offset )
        sd <- exp( beta[nb] )
        cbind( dnorm( y , mu , sd ) , pnorm( y , mu , sd ) )
    }
    ans$sim <- function( u, y, x, offset, beta) {
        nb <- length(beta)
        mu <- fm$linkinv( x %*% beta[-nb] + offset )
        sd <- exp( beta[nb] )
        qnorm( u , mu , sd )
    }
    ans$type <- "numeric"
    class(ans) <- c( "marginal.gcmr")
    ans
}

# Binomial
# y is [#success,#failure]
bn.marg <- function(link = "logit") {
    fm <- binomial( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x, offset) {
        res <- coef(glm.fit( x , y, offset=offset, family=fm ))
        names(res) <- dimnames(as.matrix(x))[[2L]]
        res
    }
    ans$npar <- function(x) NCOL(x)
    ans$dp <- function(y, x, offset, beta) {
        mu <- fm$linkinv( x %*% beta + offset )
        cbind(dbinom( y[,1] , y[,1]+y[,2], mu ) ,
              pbinom( y[,1] , y[,1]+y[,2] , mu ) )
    }
    ans$sim <- function(u, y, x, offset, beta) {
        mu <- fm$linkinv( x %*% beta + offset )
        sim <- qbinom( u, y[,1]+y[,2] , mu )
        cbind( sim, y[,1]+y[,2]-sim )
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.gcmr")
    ans
}

# Poisson
ps.marg <- function(link = "log") {
    fm <- poisson( substitute( link ) )
    ans <- list()
    ans$start <- function(y, x, offset) {
        res <- coef(glm.fit( x , y, offset=offset, family=fm ))
        names(res) <- dimnames(as.matrix(x))[[2L]]
        res
    }
    ans$npar <- function(x) NCOL(x)
    ans$dp <- function(y, x, offset, beta) {
        mu <- fm$linkinv( x %*% beta + offset )
        cbind( dpois( y , mu ) , ppois( y , mu ) )
    }
    ans$sim <- function(u, y, x, offset, beta) {
        mu <- fm$linkinv( x %*% beta + offset )
        qpois( u , mu )
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.gcmr")
    ans
}


# Negative binomial
# var(y) = E(y) + k*E(y)^2 (k>0)
nb.marg <- function(link = "log" ) {
    fm <- poisson( substitute( link ) )
    h <- 2
    h2 <- 2 - h
    eps <- .Machine$double.eps
    ans <- list()
    ans$start <- function(y, x, offset) {
        m <- glm.fit( x , y, offset=offset, family=fm )
        mu <- fitted(m)
        res <- c( coef( m ) , max(0.01, mean(((y-mu)^2-mu)/mu^h)))
        names(res) <- c(dimnames(as.matrix(x))[[2L]],"k")
        res
    }
    ans$npar <- function(x) NCOL(x)+1
    ans$dp <- function(y, x, offset, beta) {
        nb <- length(beta)
        if (beta[nb]<=0) rep(NA,NROW(y))
        mu <- fm$linkinv( x %*% beta[-nb] + offset )
        size <- pmax(eps, mu^h2 / beta[nb])
        cbind( dnbinom( y , mu=mu , size=size) , pnbinom( y , mu=mu , size=size) )
    }
    ans$sim <- function( u, y, x, offset, beta) {
        nb <- length(beta)
        if (beta[nb]<=0) rep(NA,NROW(y))
        mu <- fm$linkinv( x %*% beta[-nb] + offset )
        size <- mu^h2 / beta[nb]
        qnbinom( u , mu=mu, size=size)
    }
    ans$type <- "integer"
    class(ans) <- c( "marginal.gcmr")
    ans
}

## Skew-normal marginal
sn.marg <- function( link="identity" ) {
    if ( !require(sn) ) stop("This requires package sn")  
    fm <- gaussian( substitute( link ) ) ##  ;-)
    eps <- sqrt(.Machine$double.eps)
    ans <- list()
    ans$start <- function(y, x, offset) {
        ## "centered" sn parameterization
        ## plus further parameterization for avoiding parameters constrains
        ## Starting point is taken from sn.ml
        qrX <- qr(x)
        r <- qr.resid(qrX, y-offset)
        s <- sqrt(sum(r^2)/length(y))
        gamma1 <- sum(r^3)/(length(y) * s^3)
        if (abs(gamma1) > 0.99 ) gamma1 <- sign(gamma1) * 0.99
        ans <- c(qr.coef(qrX, y-offset), log(s), log((gamma1+0.995)/(0.995-gamma1)) )
        names(ans) <- c(dimnames(as.matrix(x))[[2L]],"log.scale", "logit.shape")
        ans
    }
    ans$npar <- function(x) NCOL(x)+2
    ans$dp <- function(y, x, offset, beta ) {
        p <- NCOL(x)
        n <- length(y)
        ## back transformation
        beta[p+1] <- exp(beta[p+1])
        beta[p+2] <- exp(beta[p+2])
        beta[p+2] <- 0.995*(beta[p+2]-1)/(beta[p+2]+1)
        ## now go back to "direct" parameterization
        beta <- cp.to.dp( beta ) 
        mu <- fm$linkinv( x %*% beta[1:p] + offset )
        scale <- beta[p+1]
        alpha <- beta[p+2] 
        cbind(dsn( y, mu, scale, alpha ) ,
              psn( y, mu, scale, alpha) )
    }
    ans$sim <- function( u, y, x, offset, beta){
        p <- NCOL(x)
        ## back transformation
        beta[p+1] <- exp(beta[p+1])
        beta[p+2] <- exp(beta[p+2])
        beta[p+2] <- 0.995*(beta[p+2]-1)/(beta[p+2]+1)
        ## now go back to "direct" parameterization
        beta <- cp.to.dp( beta ) 
        mu <- fm$linkinv( x %*% beta[1:p] + offset )
        scale <- beta[p+1]
        alpha <- beta[p+2] 
        qsn( u, mu, scale, alpha )
    }
    ans$type <- "numeric"
    class(ans) <- "marginal.gcmr"
    ans
}


 
