#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP Q_alt_b_C(SEXP Rpars, SEXP RY, SEXP RA, SEXP Rb, SEXP RSsam, SEXP Rcons, SEXP RNsim, SEXP Rn)
{
    /* Digest the datastructures (SEXPs) from R */
    int Nsim, n;
    double *pars, *Y, *A, *b, *Ssam, *cons;
    
    pars = REAL(coerceVector(Rpars, REALSXP));
    Y = REAL(coerceVector(RY, REALSXP));
    A = REAL(coerceVector(RA, REALSXP));
    b = REAL(coerceVector(Rb, REALSXP));
    Ssam = REAL(coerceVector(RSsam, REALSXP));
    cons = REAL(coerceVector(Rcons,REALSXP));
    Nsim = INTEGER(coerceVector(RNsim,INTSXP))[0];
    n = INTEGER(coerceVector(Rn,INTSXP))[0];

    /* Create SEXP to hold the answer */
    SEXP Rtemp;
    double *temp;

    PROTECT(Rtemp = allocVector(REALSXP,1));
    temp = REAL(Rtemp);

    /* main */
    int i, j;

    temp[0] = 0;
    for (i=0; i<Nsim; i++) {
	for (j=0; j<n; j++) {
	    temp[0] = temp[0]-A[j]*Ssam[i+j*Nsim]-b[j]+Y[j]*log(A[j]*Ssam[i+j*Nsim]+b[j]);
	}
    }
    temp[0] = -temp[0]/Nsim;

    /* Create SEXP for the output list */
    SEXP Rout;
    PROTECT(Rout = allocVector(VECSXP,1));
    SET_VECTOR_ELT(Rout,0,Rtemp);

    UNPROTECT(2);
    return(Rout);
}

/*
temp <- 0
for (i in 1:Nsim){
    for (j in 1:n){
        temp <- temp - A[j]*Ssam[i+j*Nsim] - b[j]*Y[j]*log(A[j]*SSam[i+j*Nsim] + b[j])
    }
}
temp = -temp/Nsim
*/

