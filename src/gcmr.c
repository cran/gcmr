#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <float.h>

void csampler( int *ipar, double *chol, double *limits, double *llik) {
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
	    for ( r=0, xr=x, mw=s2=sup=0 ; r < m ; r++, xr += n ) {
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




void bsolver( int *ipar, double *chol, double *limits, double *llik) {
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
    PutRNGstate();
}


