####

"invpBpareto" <- function(x,taus,betas,useC=TRUE){
    if (sum(betas==0)){
	ind <- (betas==0)
	taus <- taus[-ind]
	betas <- betas[-ind]
    }
    ind <- !duplicated(taus)
    taus <- taus[ind]
    betas <- betas[ind]

    betas <- betas[order(taus)]
    taus <- sort(taus)

    if (useC){
    	res <- .Call("invpBpareto_C",x,taus,betas)[[1]]
    } else {
		res <- invpBpareto_R(x,taus,betas)
    }
    return(res)
}

"invpBpareto_R" <- function(x,taus,betas)
{
	n <- length(x)
	B <- length(taus)
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
	lims[B+1] <- 1
	for (i in 1:n){
    	for (j in 1:B){
        	if ((x[i]>=lims[j])&&(x[i]<lims[j+1])){
            	res[i] <- taus[j]*((cons[j]/(1-x[i]))^(1/betas[j]))
            	break
        	}
    	}
	}
	return(res)
}

####

"rBpareto" <- function(n,taus,betas,useC=TRUE){
    return(invpBpareto(x=runif(n),taus=taus,betas=betas,useC=useC))
}

####

"optim.betas" <- function(taus,S,useC=TRUE){
    temp <- order(taus)
    taus <- sort(taus)
    taus1 <- c(taus,max(S)+1)
    if (useC){
    	out <- .Call("optim_betas_C",taus1,S)
    	betas <- out[[1]]
    	mks <- out[[2]]
    } else {
    	out <- optim_betas_R(taus1,S)
    	betas <- out[[1]]
    	mks <- out[[2]]
    }
    ##
    betas[betas==0] <- 1
    return(list(betas=betas[temp],mk=mks[temp]))
}

"optim_betas_R" <- function(taus1,S)
{
	B = length(taus1)-1
	n = length(S)
	betas <- rep(NA,B)
	mk <- rep(0,B)
	logS <- rep(0,B)
	for (i in 1:n){
    	for (j in 1:B){
        	if ((S[i]>=taus1[j]) && (S[i]<taus1[j+1])){
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
}

####

"sample.post.S.b" <- function(Y,A,taus,betas,b,S.old,Nsim,display=T,useC=TRUE){
    if(display){
		cat("Sample from MH...")
    }
    n <- length(Y)
    tempS <- rBpareto(n=Nsim*n,taus=taus,betas=betas,useC=useC)
    if (useC){
    	out <- .Call("sample_post_S_b_C",Y,A,taus,betas,b,S.old,Nsim,tempS)[[1]]
    } else {
    	tempS <- matrix(tempS,nrow=Nsim,ncol=n)
    	out <- sample_post_S_b_R(Y,A,taus,betas,b,S.old,Nsim,tempS)
    }
    if (display){
		cat("done\n")
    }
    return(out)
}

#PROTECT(Rsam=allocMatrix(REALSXP,Nsim,n));
 #    for (i=0; i<n; i++) {
 # sam[i*Nsim] = Sold[i];
 #    }

 #    GetRNGstate();

 #    for (j=1; j<Nsim; j++) {
	# for (i=0; i<n; i++) {
	#     log_accept = dpois(Y[i],A[i]*S[j+i*Nsim]+b[i],1)-dpois(Y[i],A[i]*sam[j-1+i*Nsim]+b[i],1);
	#     if (log(runif(0.0,1.0))<=log_accept){
	# 	sam[j+i*Nsim] = S[j+i*Nsim];
	#     } else {
	# 	sam[j+i*Nsim] = sam[j-1+i*Nsim];
	#     }
	# }
 #    }

"sample_post_S_b_R" <- function(Y,A,taus,betas,b,S.old,Nsim,S,display=F)
{
	n <- length(Y)
	sam <- matrix(NA,nrow=Nsim,ncol=n)
	sam[1,] <- S.old
	for (j in 2:Nsim){
        log_accept <- dpois(Y,A*S[j,]+b,log=TRUE)-dpois(Y,A*sam[j-1,]+b,log=TRUE)
        log_U <- log(runif(n))
        sam[j,] <- ifelse(log_U<=log_accept,S[j,],sam[j-1,])
        if (display){
        	cat(paste0("Number of acceptances = ",sum(log_U<=log_accept),"\n"))
        }
	}
	return(sam)
}

####

"sample.post.U.b" <- function(Y,A,taus,betas,b,U.old,Nsim,display=T,useC=TRUE){
    if(display){
	cat("Sample from MH...")
    }
    S.old <- invpBpareto(U.old,taus,betas)
    tempU <- runif(Nsim*length(Y))
    tempS <- invpBpareto(tempU,taus,betas)
    if (useC){
    	out <- .Call("sample_post_U_b_C",Y,A,taus,betas,b,U.old,S.old,Nsim,tempS,tempU)[[1]]
    } else {
    	out <- sample_post_U_b_R(Y,A,taus,betas,b,U.old,S.old,Nsim,tempS,tempU)
    }
    if (display){
		cat("done\n")
    }
    return(out)
}

"sample_post_U_b_R" <- function(Y,A,taus,betas,b,Uold,Sold,Nsim,S,U)
{
	n <- length(Y)
	sam <- matrix(NA,Nsim,n)
	SSam <- matrix (NA,Nsim,n)
	sam[1,] <- Uold
	Ssam[1,] <- Sold
	for (j in 2:Nsim){
		log_accept <- dpois(Y,A+S[j,]+b,log=TRUE) - dpois(Y,A+Ssam[j-1,]+b,log=TRUE)
		log_U <- log(runif(n))
		sam[j,] <- ifelse(log_U<log_accept,U[j,],sam[j-1,])
		Ssam[j,] <- ifelse(log_U<log_accept,S[j,],Ssam[j-1,])
	}
	return(list(sam,Ssam))
}

 #    for (j=1; j<Nsim; j++) {
	# for (i=0; i<n; i++) {
	#     accept = exp(dpois(Y[i],A[i]*S[j+i*Nsim]+b[i],1)-dpois(Y[i],A[i]*Ssam[j-1+i*Nsim]+b[i],1));
	#     if (unif_rand()<=accept){
	# 	sam[j+i*Nsim] = U[j+i*Nsim];
	# 	Ssam[j+i*Nsim] = S[j+i*Nsim];
	#     } else {
	# 	sam[j+i*Nsim] = sam[j-1+i*Nsim];
	# 	Ssam[j+i*Nsim] = Ssam[j-1+i*Nsim];
	#     }
	# }
 #    }

####

"Q.alt.b" <- function(pars,Y,A,b,sam,const,useC=TRUE){
    #cat(pars,const,"\n")
    B <- length(pars)/2
    pars[1:B] <- pars[1:B]/const
    S.sam <- matrix(invpBpareto(as.vector(sam),pars[1:B],pars[(B+1):(2*B)]),nrow=nrow(sam))
    if (useC){
    	out <- .Call("Q_alt_b_C",pars,Y,A,b,S.sam,const,nrow(S.sam),ncol(S.sam))
    	temp <- out[[1]]
    } else {
    	temp <- Q_alt_b_R(pars,Y,A,b,S.sam,const,nrow(S.sam),ncol(S.sam))
    }
    if (is.na(temp)||is.infinite(temp)){
	return(9e300)
    } else {
	return(temp)
    }
}

"Q_alt_b_R" <- function(pars,Y,A,b,Ssam,const,Nsim,n)
{
	temp <- 0
	for (i in 1:Nsim){
    	for (j in 1:n){
        	temp <- temp - A[j]*Ssam[i+j*Nsim] - b[j]*Y[j]*log(A[j]*Ssam[i+j*Nsim] + b[j])
    	}
	}
	temp = -temp/Nsim
	return(temp)
}

####

"sample.power.post.S.b" <- function(Y,A,taus,betas,b,S.old,pwr,Nsim,display=T,useC=TRUE){
    if(display){
		cat("Sample from MH...")
    }
    n <- length(Y)
    tempS <- rBpareto(n=Nsim*n,taus=taus,betas=betas,useC=useC)
    if (useC){
    	out <- .Call("sample_power_post_S_b_C",Y,A,taus,betas,b,S.old,pwr,Nsim,tempS)[[1]]
    } else {
    	tempS <- matrix(tempS,nrow=Nsim,ncol=n)
    	out <- sample_power_post_S_b_R(Y,A,taus,betas,b,S.old,pwr,Nsim,tempS)
    }
    if (display){
		cat("done\n")
    }
    return(out)
}

"sample_power_post_S_b_R" <- function(Y,A,taus,betas,b,Sold,pwr,Nsim,S)
{
    n <- length(Y)
    sam <- matrix(NA,Nsim,n)
    sam[1,] <- Sold
    for (j in 2:Nsim){
        log_accept <- pwr*(dpois(Y,A+S[j,]+b,log=TRUE)-dpois(Y,A+sam[j-1,]+b,log=TRUE))
        log_U <- log(runif(n))
        sam[j,] <- ifelse(log_U<=log_accept,S[j,],sam[j-1,])
    }
    return(sam)
}
