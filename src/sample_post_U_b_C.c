#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP sample_post_U_b_C(SEXP RY, SEXP RA, SEXP Rtaus, SEXP Rbetas, SEXP Rb, SEXP RUold, SEXP RSold, SEXP RNsim, SEXP RS, SEXP RU)
{
    /* Digest the datastructures (SEXPs) from R */
    int Nsim, n;
    double *Y, *A, *taus, *betas, *b, *Uold, *Sold, *U, *S;

    Y = REAL(coerceVector(RY,REALSXP));
    A = REAL(coerceVector(RA,REALSXP));
    taus = REAL(coerceVector(Rtaus,REALSXP));
    betas = REAL(coerceVector(Rbetas,REALSXP));
    b = REAL(coerceVector(Rb,REALSXP));
    Uold = REAL(coerceVector(RUold,REALSXP));
    Sold = REAL(coerceVector(RSold,REALSXP));
    U = REAL(coerceVector(RU,REALSXP));
    S = REAL(coerceVector(RS,REALSXP));
    Nsim = INTEGER(coerceVector(RNsim,INTSXP))[0];
    n = length(RY);

    /* Create SEXP to hold the answer */
    SEXP Rsam, RSsam;
    double *sam, *Ssam;
    
    PROTECT(Rsam=allocMatrix(REALSXP,Nsim,n));
    PROTECT(RSsam=allocMatrix(REALSXP,Nsim,n));
    sam = REAL(Rsam);
    Ssam = REAL(RSsam);

    /* main */
    int i, j;
    double log_accept;
    
    for (i=0; i<n; i++) {
	sam[i*Nsim] = Uold[i];
	Ssam[i*Nsim] = Sold[i];
    }

    GetRNGstate();

    for (j=1; j<Nsim; j++) {
	for (i=0; i<n; i++) {
	    log_accept = dpois(Y[i],A[i]*S[j+i*Nsim]+b[i],1)-dpois(Y[i],A[i]*Ssam[j-1+i*Nsim]+b[i],1);
	    if (log(runif(0.0,1.0))<=log_accept){
		sam[j+i*Nsim] = U[j+i*Nsim];
		Ssam[j+i*Nsim] = S[j+i*Nsim];
	    } else {
		sam[j+i*Nsim] = sam[j-1+i*Nsim];
		Ssam[j+i*Nsim] = Ssam[j-1+i*Nsim];
	    }
	}
    }

    PutRNGstate();

    /* Create SEXP for the output list */
    SEXP Rout;
    PROTECT(Rout = allocVector(VECSXP,2));
    SET_VECTOR_ELT(Rout,0,Rsam);
    SET_VECTOR_ELT(Rout,1,RSsam);

    UNPROTECT(3);
    return(Rout);
}

/*
sam <- matrix(NA,Nsim,n)
SSam <- matrix (NA,Nsim,n)
sam[1,] <- Uold
Ssam[1,] <- Sold
for (j in 2:Nsim){
    for (i in 1:n){
        log_accept <- dpois(Y[i],A[i]+S[j,i]+b[i],log=TRUE) - dpois(Y[i],A[i]+Ssam[j-1,i]+b[i],log=TRUE)
        if (log(runif(1))<log_accept){
            sam[j,i] <- U[j,i]
            Ssam[j,i] <- S[j,i]
        } else {
            sam[j,i] <- sam[j-1,i]
            Ssam[j,i] <- Ssam[j-1,i]
        }
    }
}
return(list(sam,Ssam))
*/
