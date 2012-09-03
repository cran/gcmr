#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <float.h>

/* Forward declarations */
static void csampler( int *ipar, double *chol, double *limits, double *llik);
static void bsolver( int *ipar, double *chol, double *limits, double *llik);

/* 
   R interface.

   Input:   
   - ipar=c(is response discrete?,should I return the likelihood?,
            number of observations, number of Monte Carlo replications,
            number of clusters, length of cluster 1, lenght of cluster 2,...)
   - chol=Cholesky factor of the correlation matrix.
          If the number of clusters is greater than 1, it can be
          the Cholesky factor of the correlation matrix within the first cluster,
	  followed by the Cholesky factor of the correlation matrix within the second,
	  and so on.
   - dp = a nx2 matrix containing the value of the marginal density and cdf

   Output:
   - if (should I return the likelihood != 0) a vector with an
     approximation of log p(y_i|y_{i-1},...,y_1).
   - otherwise a nx2 matrix containing what is needed
     for the computation of the conditional residuals. 
     In particular, if the response is discrete, 
     the matrix contains an approximation of
     [first column] p(y_i|y_{i-1},...,y_1)
     [second column] p(Y_i<=y_i|Y_{i-1}=y_{i-1},...,Y_{1}=y_1),
     i.e., a conditional version of dp.
     If the response is continuous, the second column contains
     the conditional residuals.
*/
static SEXP gcmrcomp(SEXP ipar, SEXP chol, SEXP dp) {
    SEXP llik;
    /* coerce the arguments to the right type */
    PROTECT(ipar = coerceVector(ipar, INTSXP));
    PROTECT(chol = coerceVector(chol, REALSXP));
    /* Since dp is modified in csampler/bsolver we always make a fresh copy */
    PROTECT(dp = (TYPEOF(dp)==REALSXP) ? duplicate(dp) : coerceVector(dp,REALSXP));
    /* allocate the space for the log likelihood */
    PROTECT(llik = allocVector(REALSXP, INTEGER(ipar)[2]));
    /* Ok, we can do the real thing... */
    if (INTEGER(ipar)[0]) {
	csampler(INTEGER(ipar)+2,REAL(chol),REAL(dp),REAL(llik));
    } else {
	bsolver(INTEGER(ipar)+2,REAL(chol),REAL(dp),REAL(llik));
    }
    /* and return */
    UNPROTECT(4);
    return((INTEGER(ipar)[1]) ? llik : dp);
} 


/* Registration */
static const R_CallMethodDef CallMethods[] = {
    {"gcmrcomp", (DL_FUNC) &gcmrcomp, 3},
    {NULL,NULL,0}
};

void R_init_gcmr(DllInfo *info) {
    R_registerRoutines(info, NULL, CallMethods, NULL, NULL);
}


/* 
   Likelihood and conditional residuals computation in the discrete case:
   Conditional (GHK) sampler 
*/
static void csampler( int *ipar, double *chol, double *limits, double *llik) {
    int g, ig, ng, r , i, j , ij , n = ipar[0] , m = ipar[1] , nstrata = ipar[2],
	*lstrata = ipar + 3;
    double *a = limits, *b = limits+n, ZERO=0, ONE=1, yi, cii , sup,
	lk , plow, pup, z, mw, s2, mwold, biasold, EPS = sqrt(DOUBLE_EPS), EPS1=1-EPS,
	*w = (double *) R_alloc(n*m+m,sizeof(double)), *x=w+m, *xr;
    GetRNGstate();
    for ( g=0, i=0 ; g < nstrata ; g++) {
	ng = lstrata[g] ;
	mwold = 1 ;
	biasold = 0 ;
	for ( r=0 ; r<m ; r++) w[r] = ONE;
	for ( ig=0 ; ig<ng ; ig++, i++) {
	    a[i] = qnorm(fmax2(EPS,fmin2(EPS1,b[i]-a[i])), ZERO, ONE, 1, 0);
	    b[i] = qnorm(fmax2(EPS,fmin2(EPS1,b[i])), ZERO, ONE, 1, 0);
	    mw = s2 = sup = 0.0 ;
	    for ( r=0 ; r < m ; r++ ) {
		xr = x + r*n ;
		for ( j=0 , yi=0 , ij = ig*ng ; j < ig ; j++ , ij++) yi += chol[ij]*xr[j] ;
		cii = chol[ig+ig*ng] ;
		plow = pnorm(a[i] , yi , cii, 1 , 0 ) ;
		pup = pnorm(b[i] , yi , cii, 1 , 0 ) ;
		z = fmax2(EPS,fmin2(EPS1,plow+unif_rand()*(pup-plow)));
		z = qnorm( z , yi , cii , 1 , 0 ) ;
		xr[ig] = ( z - yi ) / cii ;
		w[r] /= mwold ;
		sup += pup*w[r] ;
		lk = fmax2(EPS,fmin2(EPS1,pup-plow));
		w[r] *= lk ;
		mw += w[r] ;
		s2 += w[r]*w[r];
	    }
	    mw /= m ;
	    s2 = ((s2/m)-mw*mw)/(2*(m-1)*mw*mw) ;
	    llik[i] = log(mw) + s2 - biasold ;
	    mwold = mw ;
	    biasold = s2 ;
	    a[i] = mw ;
	    b[i] = sup / m ;
	}
	x += ng ;
	chol += ng*ng ;
    }
    PutRNGstate();
}



/* 
   Likelihood and residuals computation in the continuos case 
*/
static void bsolver( int *ipar, double *chol, double *limits, double *llik) {
    int g, ig, ng, i, j , ij , n = ipar[0] , m = ipar[1] , nstrata = ipar[2],
	*lstrata = ipar + 3;
    double *a = limits, *b = limits+n, ZERO=0, ONE=1, TWO=2, yi, cii , lk , z, 
	*es=b , EPS = sqrt(DOUBLE_EPS), EPS1=1-EPS;
    for ( g=0, i=0 ; g<nstrata ; g++) {
	ng = lstrata[g] ;
	for ( ig=0 ; ig<ng ; ig++, i++) {
	    for ( j=0 , yi=0 , ij = ig*ng ; j < ig ; j++ , ij++) yi += chol[ij]*es[j] ;
	    z = qnorm(fmax2(EPS,fmin2(EPS1,b[i])), ZERO, ONE, 1, 0);
	    cii = chol[ig+ig*ng] ;
	    b[i] = (z - yi) / cii ;
	    llik[i] = (z*z-b[i]*b[i])/TWO+log(a[i]/cii) ;
	}
	es += ng ;
	chol += ng*ng ;
    }
}

