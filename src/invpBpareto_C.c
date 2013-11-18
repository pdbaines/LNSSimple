#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP invpBpareto_C(SEXP Rx, SEXP Rtaus, SEXP Rbetas)
{
    /* Digest the datastructures (SEXPs) from R */
    int B, n;
    double *x, *taus, *betas;

    x = REAL(coerceVector(Rx,REALSXP));
    taus = REAL(coerceVector(Rtaus,REALSXP));
    betas = REAL(coerceVector(Rbetas,REALSXP));
    B = length(Rtaus);
    n = length(Rx);

    /* Create SEXP to hold the answer */
    SEXP Rres;
    double *res;

    PROTECT(Rres=allocVector(REALSXP,n));
    res = REAL(Rres);

    /* main */
    int i, j;
    double cons[B], lims[B+1];

    cons[0] = 1;
    lims[0] = 0;
    for (j=1; j<B; j++) {
	cons[j] = cons[j-1]*pow(taus[j-1]/taus[j],betas[j-1]);
	lims[j] = 1-cons[j];
    }
    lims[B] = 1;
    
    for (i=0; i<n; i++) {
	for (j=0; j<B; j++) {
	    if ((x[i]>=lims[j])&&(x[i]<lims[j+1])) {
		res[i] = taus[j]*pow(cons[j]/(1-x[i]),1/betas[j]);
		break;
	    }
	}
    }

    /* Create SEXP for the output list */
    SEXP Rout;
    PROTECT(Rout=allocVector(VECSXP,1));
    SET_VECTOR_ELT(Rout,0,Rres);

    UNPROTECT(2);
    return(Rout);
}

/*
res <- rep(NA,n)
cons <- rep(NA,B)
lims <- rep(NA,B+1)
cons[1] <- 1
lims[1] <- 0
if (B>1){
    for (j in 2:B){
        cons[j] <- cons[j-1]*((taus[j-1]/taus[j])^(betas[j-1]))
        lims[j] <- 1-cons[j]
    }
}
lims[B] <- 1
for (i in 1:n){
    for (j in 1:B){
        if ((x[i]>=lims[j])&&(x[i]<lims[j+1])){
            res[i] <- taus[j]*((cons[j]/(1-x[i]))^(1/betas[j]))
            break
        }
    }
}
return(res)
*/
