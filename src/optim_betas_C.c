#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP optim_betas_C(SEXP Rtaus1, SEXP RS)
{
    /* Digest the datastructures (SEXPs) from R */
    int B, n;
    double *taus1, *S;

    taus1 = REAL(coerceVector(Rtaus1,REALSXP));
    S = REAL(coerceVector(RS,REALSXP));
    B = length(Rtaus1)-1;
    n = length(RS);

    /* Create SEXP to hold the answer */
    SEXP Rbetas, Rmk;
    int *mk;
    double *betas;

    PROTECT(Rbetas=allocVector(REALSXP,B));
    PROTECT(Rmk=allocVector(INTSXP,B));
    betas = REAL(Rbetas);
    mk = INTEGER(Rmk);

    /* main */
    int i, j, nk, nk1;
    double logS[B];

    for (j=0; j<B; j++) {
	logS[j] = 0;
	mk[j] = 0;
    }

    for (i=0; i<n; i++) {
	for (j=0; j<B; j++) {
	    if ((S[i]>=taus1[j])&&(S[i]<taus1[j+1])) {
		mk[j] = mk[j]+1;
		logS[j] = logS[j]+log(S[i]);
		break;
	    }
	}
    }

    nk = n;
    for (j=0; j<B; j++) {
	nk1 = nk-mk[j];
	betas[j] = mk[j]/(logS[j]+nk1*log(taus1[j+1])-nk*log(taus1[j]));
	nk = nk1;
    }

    /* Create SEXP for the output list */
    SEXP Rout;
    PROTECT(Rout = allocVector(VECSXP,2));
    SET_VECTOR_ELT(Rout,0,Rbetas);
    SET_VECTOR_ELT(Rout,1,Rmk);

    UNPROTECT(3);
    return(Rout);
}

/*
B = length(taus1)-1
n = length(S)
betas <- rep(NA,B)
mk <- rep(0,B)
logS <- rep(0,B)
for (i in 1:n){
    for (j in 1:B){
        if ((S[i]>=taus1[j])&&(S[i]<taus1[j+1])){
            mk[j] = mk[j]+1
            logS[j] = logS[j]+log(S[i])
            break
        }
    }
}
nk <- n
for (j in 1:B){
    nk1 <- nk-mk[j]
    betas[j] <- mk[j]/(logS[j]+nk1*log(taus1[j+1])-nk*log(taus1[j]))
    nk <- nk1
}
return(list(betas,mk))
*/
