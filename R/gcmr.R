# Return a function which computes an approximation of log(c(p(y[1]),p(y[2]|y[1]),...))
llik <- function(m, irep, magic=NA) {
    n <- m$n
    not.na <- m$not.na
    y <- m$y[not.na,]
    x <- m$x[not.na,,drop=FALSE]
    offset <- m$offset[not.na,]
    is.int <- m$marginal$type == "integer"
    ind.lik <- !is.null(m$cormat$independent) 
    dp <- m$marginal$dp
    ibeta <- m$ibeta
    igamma <- m$igamma
    vchol <- m$cormat$chol
    if (missing(irep)) irep <- max(m$options$nrep)
    seed <- m$options$seed
    id <- if(is.null(m$cormat$id)) rep(1,n) else m$cormat$id
    lstrata <- rle(id)$length
    ipar <- c(if(is.int) 1 else 0,1,n,irep,length(lstrata),lstrata)
    magic <- rep(magic,m$n)
    theta <- m$fixed
    ifree <- is.na(theta)
    cache <- new.env()
    function( theta.free ) {
        theta[ifree] <- theta.free
        beta <- theta[ibeta]
        if (!identical(cache$beta,beta)) {
            dp <- dp(y,x,offset,beta)
            assign("beta",beta,envir=cache)
            assign("dp",dp,envir=cache)
        } else {
            dp <- get("dp",envir=cache)
        }
        if ( is.null(dp) ) return( magic )
        if ( ind.lik ) return( log(dp[,1]) )
        gamma <- theta[igamma]
        if (!identical(cache$gamma,gamma)) {
            q <- vchol( gamma , not.na )
            assign("gamma",gamma,envir=cache)
            assign("q",q,envir=cache)
        } else {
            q <- get("q",envir=cache)
        }
        if ( is.null(q) ) return ( magic )
        if (is.int) set.seed(seed)
        lk <- .Call(gcmrcomp, ipar, unlist(q), dp)
        if ( all(is.finite(lk)) ) lk else magic
    }    
} 

# Workhorse which is called from gcmr.fit and profile.gcmr
truefit <- function(x) {
    # saving/restoring the random seed
    if ( x$marginal$type == "integer" ) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
        }
    }
    ifree <- is.na( x$fixed )
    theta <- x$fixed
    theta[ifree] <- x$estimate[ifree]
    low <- x$lower[ifree]
    up <- x$upper[ifree]
    big <- -sqrt(.Machine$double.xmax)
    nrep <- if(x$marginal$type=="numeric") 1 else x$options$nrep
    for ( i in 1:length(nrep) ) {
        log.lik <- llik(x,nrep[i],big)
        if( i!=length(nrep) ) # no warnings until last optimization
          ans <- suppressWarnings( x$options$opt(theta[ifree] , log.lik , low, up) )
        else
          ans <- x$options$opt( theta[ifree] , log.lik , low, up)
        theta[ifree] <- ans$estimate
    }
    names(theta) <- names(x$estimate)
    x$estimate <- theta
    x$maximum <- ans$maximum
    x
}

# Add/replace an estimate of jacobian and hessian to an gcmr object
# Hessian is approximated by finite differences
# in a rotated space in which the true hessian is near
# the identity matrix. The rotation is based on
# jacobian crossprod
jhess <- function(x, options=x$options, only.jac = FALSE) {
    if ( !inherits( x , "gcmr" ) ) stop("First argument must be a gcmr object")
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
    }
    theta <- x$estimate
    ifree <- is.na( x$fixed )
    theta.free <- theta[ifree]
    low <- x$lower[ifree]
    up <- x$upper[ifree]
    xlik <- llik(x)
    log.lik <- function(th) xlik(pmax(low,pmin(up,th)))
    eps <- .Machine$double.eps^(1/4)
    relStep <- 0.1
    maxtry <- 10
    delta <- ifelse(abs(theta.free)<1, eps, eps*theta.free)
    di <- function(i,delta) {
        x1 <- x2 <- theta.free
        x1[i] <- x1[i] - delta[i]
        x2[i] <- x2[i] + delta[i]
        (log.lik(x2)-log.lik(x1))/(2*delta[i])
    }
    while (1) {
        x$jac <- sapply(seq_along(theta.free),di,delta)
        if( all(is.finite(x$jac)) ) break
        delta <- delta/2
        maxtry <- maxtry - 1
        if (maxtry<0) stop("impossible to compute a finite jacobian")
    }
    if (!only.jac) {
        a <- svd(x$jac)
        a$d <- pmax(a$d,sqrt(.Machine$double.eps)*a$d[1])
        x$hessian <- nlme:::fdHess(rep(0,length(theta.free)),
                                   function(tx) sum(log.lik(theta.free+a$v%*%(tx/a$d))),
                                   minAbsPar=1,.relStep=relStep)$Hessian
        x$hessian = (x$hessian+t(x$hessian))/2
        x$hessian <- a$v%*%(outer(a$d,a$d)*(x$hessian))%*%t(a$v)
    }
    x$options$no.se <- FALSE
    x
}

gradcheck <- function(x, options=x$options) {
    j <- jhess(x,options,only.jac=TRUE)$jac
    g <- colSums(j)
    sum(g*solve(crossprod(j),g))
}

gcmr.opt <- function(start,loglik, lower, upper) {
    ans <- nlminb(start,function(x) -sum(loglik(x)), lower=lower, upper=upper)
    if(ans$convergence) warning(paste("nlminb exits with code",ans$convergence))
    list(estimate=ans$par,maximum=ans$objective)
}

gcmr.options <- function(seed=round(runif(1,1,100000)), nrep=c(100,1000),
                         no.se=FALSE, opt=gcmr.opt) {
    list(seed=seed,nrep=nrep,no.se=no.se,opt=opt)
}

gcmr <- function(formula, data, subset, offset, contrasts=NULL, marginal, cormat,
                 start, fixed, options=gcmr.options()){
    call <- match.call()
    if( !missing(data) ) {
        marginal <- eval(call$marginal,data,parent.frame())
        cormat <- eval(call$cormat,data,parent.frame())
    }
    if( is.function( marginal ) )
     marginal <- marginal()
    if ( !inherits( marginal , "marginal.gcmr" ) ) stop("Unknown marginal")
    if( is.function( cormat ) )
      cormat <- cormat()
    if ( !inherits( cormat , "cormat.gcmr" ) ) stop("Unknown cormat")
    if (missing(data)) 
        data <- environment(formula)
    ## next lines partially "inherited" from function glm
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- na.pass
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
    if (is.null(offset)) 
        offset <- rep.int(0, NROW(Y))
    fit <- gcmr.fit(X, Y, offset, marginal, cormat, start, fixed, options)
    fit$call <- call
    fit
}

gcmr.fit <- function(x=rep(1,NROW(y)), y, offset=rep(0,NROW(y)), marginal, cormat, start, fixed, options=gcmr.options()) {
   if( is.function( marginal ) )
     marginal <- marginal()
   if( is.function( cormat ) )
     cormat <- cormat()
    ## arguments check
    if ( !inherits( marginal , "marginal.gcmr" ) ) stop("Unknown marginal")
    if ( !inherits( cormat , "cormat.gcmr" ) ) stop("Unknown cormat")
    ## gcmr object
    nb <- marginal$npar(x)
    ng <- cormat$npar
    not.na <- apply(cbind(y,x,offset),1,function(z) !any(is.na(z)))
    m <- structure(list(y=as.matrix(y), x=as.matrix(x), offset=as.matrix(offset),
                        not.na=not.na, n=sum(not.na), marginal=marginal, cormat=cormat,
                        ibeta=1:nb, igamma=if (ng) (nb+1):(nb+ng) else NULL,
                        call=match.call()), class="gcmr")
    if ( missing(start) ) {
        lambda <- marginal$start(m$y[not.na,],m$x[not.na,,drop=FALSE],m$offset[not.na,])
        tau <- cormat$start()
        m$estimate <- c(lambda,tau)
        ll <- attr(lambda,"lower")
        lt <- attr(tau,"lower")
        m$lower <- c(if(is.null(ll)) rep(-Inf,length(lambda)) else ll,
                     if(is.null(lt)) rep(-Inf,length(tau)) else lt)
        ll <- attr(lambda,"upper")
        lt <- attr(tau,"upper")
        m$upper <- c(if(is.null(ll)) rep(Inf,length(lambda)) else ll,
                     if(is.null(lt)) rep(Inf,length(tau)) else lt)
    } else {
        m$estimate <- start
        ll <- attr(start,"lower")
        m$lower <- if(is.null(ll)) rep(-Inf,length(start)) else ll
        ll <- attr(start,"upper")
        m$upper <- if(is.null(ll)) rep(Inf,length(start)) else ll
    }
    if (length(m$estimate) != nb+ng) stop("mismatch in the number of initial parameters")
    if (length(m$lower) != nb+ng) stop("length of lower different from the number of the parameters")
    if (length(m$upper) != nb+ng) stop("length of upper different from the number of the parameters")
   if ( missing(fixed) ) {
        m$fixed <- rep( NA , length(m$estimate) )
    } else {
        if (length(fixed) != length(m$estimate) ) stop("fixed has a wrong length")
        m$fixed <- fixed
    }
    m$options <- do.call(gcmr.options,options)
    # compute estimate
    m <- truefit(m)
    # and s.e. and return
    if (m$options$no.se) m else jhess(m)
}

coef.gcmr <- function(object,...) object$estimate

logLik.gcmr <- function(object,...) {
    ans <- -object$maximum
    attr(ans,"df") <- sum( is.na(object$fixed) )
    class(ans) <- "logLik"
    ans
}

pinv <- function(h,tol) {
    h <- svd(h)
    idx <- h$d > sqrt(.Machine$double.eps)*h$d[1]
    ans <- h$u[,idx,drop=FALSE]%*%( (1/h$d[idx])*t(h$u[,idx,drop=FALSE]))
    attr(ans,"df") <- sum(idx)
    ans
}

estfun.gcmr <- function(x,...) x$jac
bread.gcmr <- function(x,...) pinv(x$hessian)*x$n

vcov.gcmr <- function(object, type=c("hessian","sandwich","vscore","cluster","hac"),...) {
    if ( !inherits( object , "gcmr" ) ) stop("First argument must be a mr object")
    type <- match.arg(type)
    if (type=="sandwich") {
        if (is.null(object$hessian) || is.null(object$jac)) object <- jhess(object)
        h <- pinv(object$hessian)
        v <- h %*% crossprod(object$jac) %*% h
    } else if (type=="hessian")  {
        if (is.null(object$hessian) || is.null(object$jac)) object <- jhess(object)
        v <- pinv(object$hessian)
    } else if (type=="vscore")  {
        if (is.null(object$jac)) object <- jhess(object,only.jac=TRUE)
        v <- pinv(crossprod(object$jac))
    } else if (type=="cluster") {
        if (is.null(object$hessian) || is.null(object$jac)) object <- jhess(object)
        if (is.null(object$cormat$id)) stop("no cluster found")
        h <- pinv(object$hessian)
        v <- sapply(split(1:object$n,object$cormat$id),function(i) colSums(object$jac[i,])) 
        v <- h %*% (v %*% t(v) ) %*% h
    } else {
        if (!require(sandwich)) stop("HAC requires package sandwich")
        if (is.null(object$hessian) || is.null(object$jac)) object <- jhess(object)
        v <- vcovHAC(object,...)
    }
    ans <- matrix(0,length(object$estimate),length(object$estimate))
    ifree <- is.na(object$fixed)
    ans[ifree,ifree] <- v
    colnames(ans) <- rownames(ans) <- names(object$estimate)
    ans
}

se <- function(x,...) sqrt(diag(vcov.gcmr(x,...)))

print.gcmr <- function (x, digits = max(3, getOption("digits") - 3), k = 2 ,...)
{
    cat("\nCall:", deparse(x$call, width.cutoff = 75), "", sep = "\n")
    cat("Parameters:\n")
    par <- x$estimate
    stderr <- if (x$options$no.se) NULL else se(x,...)
    xpar <- round( rbind( as.numeric(par) , s.e.=stderr ) , digits=digits )
    colnames(xpar) <- names(par)
    print.default(xpar, print.gap = 2)
    cat("\nlog likelihood = ", format(round(-x$maximum, 2)),
        ",  aic = ", format(round(AIC(x,k=k), 2)), "\n", sep = "")
    invisible(x)
}

## summary/print to be improved in future
summary.gcmr <- function(object, ...)
  print.gcmr(object)

## note: method=interval is hidden
residuals.gcmr <- function (object, type=c("conditional","marginal"),
                          method=c("random","mid"),...) {
    type <- match.arg(type)
    method <- match.arg(method)
    is.int <- object$marginal$type == "integer"
    cond <- is.null(object$cormat$independent) && (type=="conditional") 
    if (is.int) {
        ##before saving/setting seed to obtain different randomization 
        if (method=="random") u <- runif(object$n) 
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            seed.keep <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
        }
        set.seed(object$options$seed)
    }
    a <- object$marginal$dp(object$y[object$not.na,],
                            object$x[object$not.na,,drop=FALSE],
                            object$offset[object$not.na,],
                            object$estimate[object$ibeta])
    if (cond) {
        q <- object$cormat$chol(object$estimate[object$igamma],object$not.na)
        id <- if(is.null(object$cormat$id)) rep(1,object$n) else object$cormat$id
        lstrata <- rle(id)$length
        ipar <- c(if(is.int) 1 else 0, 0,
                  object$n,max(object$options$nrep),length(lstrata),lstrata)
        a <- .Call(gcmrcomp,ipar,unlist(q),a)
    }
    if (is.int) {
        if (method=="random") {
            ans <- rep(NA,length(object$not.na))
            ans[object$not.na] <- qnorm(a[,2]-u*a[,1])            
        } else if (method=="interval") {
            ans <- matrix(NA,length(object$not.na),2)
            ans[object$not.na,] <- qnorm(cbind(a[,2],a[,2]-a[,1]))
        } else {
            ans <- rep(NA,length(object$not.na))
            ans[object$not.na] <- qnorm(a[,2]-a[,1]/2)            
        }
    } else {
            ans <- rep(NA,length(object$not.na))
            ans[object$not.na] <- if (cond) a[,2] else qnorm(a[,2])            
    }
    ans
}

plotint <- function(r) {
    r <- as.ts(r)
    tm <- time(r)
    plot(NA,NA,xlim=range(tm),xlab="",ylim=range(r),ylab="pseudo residuals")
    for (i in 1:length(tm)) lines(c(tm[i],tm[i]),r[i,])
}

profile.gcmr <- function(fitted , which , low , up, npoints = 10 , display = TRUE , alpha = 0.05, ... ) {
    if ( !inherits( fitted , "gcmr" ) ) stop("first argument must be a mr object")
    if(is.null(low) || is.null(up)){
      this.se <- se(fitted)[which]
      if(missing(low))
        low <- coef(fitted)[which]-3*this.se
      if(missing(up))
        up <- coef(fitted)[which]+3*this.se
    }
    points <- seq( low , up , length = npoints )
    prof <- function(x,which,fitted) {
        fitted$fixed[which] <- x
        gcmr:::truefit(fitted)$maximum
    }
    loglik <- sapply(points,prof,which,fitted)
    points <- c( points , fitted$par[which] )
    ord <- order( points)
    points <- points[ord]
    loglik <- -c(loglik,fitted$maximum)[ord]
    if ( display ) {
        npoints <- seq( low , up , length = 200 )
        plot( npoints , splinefun( points, loglik )(npoints) , type = "l" ,
             xlab = names(fitted$estimate)[which], ylab = "log-likelihood profile")
        grid()
        abline( h = max(loglik) - qchisq( 1 - alpha , 1 )/2 , lty = "dashed" )
    }
    invisible(list(points=points,profile=loglik))
}

hausman <- function(m, method=c("one-step","full"),
                      nboot=if (method=="one-step") 1000 else 200,
                      nrep=200) {
    if ( !inherits( m , "gcmr" ) ) stop("first argument must be a gcmr object")
    method <- match.arg(method)
    m$options$nrep <- nrep
    m$options$no.se <- TRUE
    ## ind lik
    mind <- gcmr.fit(x=m$x,y=m$y,offset=m$offset,marginal=m$marg,cormat=ind.cormat(),fixed=m$fixed[m$ibeta],options=m$options)
    ## differences in the marginal parameters
    betaml <- coef(m)[m$ibeta]
    betaind <- mind$estimate
    delta <- betaind - betaml
    ## simulator
    R <- m$cormat$chol( coef(m)[m$igamma] , m$not.na )
    if (!is.list(R)) R <- list(R)
    R <- lapply(R,t)
    y <- m$y[m$not.na,]
    x <- m$x[m$not.na,,drop=FALSE]
    offset <- m$offset[m$not.na,]
    sim <- function() {
        u <- pnorm(unlist(lapply(R, function(x) x %*% rnorm(NROW(x)))))
        m$marginal$q(u,x,offset,betaml)
    }
    ## bootstrap estimate of the the variance of the difference
    vdelta <- var(if (method=="one-step") dh1(m,mind,sim,nboot) else dhf(m,mind,sim,nboot))
    ## return
    stderr <- sqrt(diag(vdelta))
    htab <- cbind(ml=betaml,ind.lik=betaind,diff=delta,s.e.=stderr,z=delta/stderr)
    rownames(htab) <- names(betaml)
    vdelta <- pinv(vdelta)
    hall <- sum(delta*(vdelta%*%delta))
    hall <- c(stat=hall,df=attr(vdelta,"df"),
              p.value=pchisq(hall,attr(vdelta,"df"),lower.tail=FALSE))
    list(overall=hall,parameters=htab,
         options=c(method=method,nboot=nboot,nrep=nrep))    
}

dh1 <- function(m, mind, sim, nboot) {
    mind$estimate <- coef(m)[m$ibeta]
    nb <- Hind <- 0
    sml <- matrix(0,nboot,length(coef(m))) 
    sind <- matrix(0,nboot,length(coef(mind))) 
    while ( nb < nboot ) {
        nb <- nb+1
        m$y[m$not.na,] <- mind$y[mind$not.na,] <- sim()
        sml[nb,] <- colSums(jhess(m, only.jac=TRUE)$jac)
        mind <- jhess(mind, only.jac=FALSE)
        sind[nb,] <- colSums(mind$jac) 
        Hind <- Hind + (mind$hessian-Hind)/nb
    }
    (sml%*%pinv(crossprod(sml/sqrt(nboot))))[,m$ibeta]-sind%*%pinv(Hind)
}

dhf <- function(m, mind, sim, nboot) {
    theta <- m$estimate
    m$options$opt <- mind$options$opt <- function(start,loglik) {
        ans <- nlminb(start,function(x) -sum(loglik(x)),
                      control=list(rel.tol=1e-4,iter.max=3*length(theta)))
        list(estimate=ans$par,maximum=ans$objective)
    }
    d <- matrix(0,nboot,length(coef(mind)))
    nb <- 0
    while (nb < nboot) {
        nb <- nb+1
        m$y[m$not.na,] <- mind$y[mind$not.na,] <- sim()
        m$estimate <- theta
        m <- truefit(m)
        mind$estimate <- mind$marginal$start(mind$y,mind$x,mind$offset)
        mind <- truefit(mind)
        d[nb,] <- coef(m)[m$ibeta]-coef(mind)
    }
    d
}
