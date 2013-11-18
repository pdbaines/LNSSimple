#######################
### initial values ####
#######################
# logN-logS plot
logNlogS <- function(S){
  S <- sort(S)
  n <- length(S)
  vals <- unique(S)

  logN <- log10(n-cumsum(tabulate(match(S, vals))))[-length(vals)]
  logS <- log10(vals)[-length(vals)]
  #rval <- approxfun(logS,logN,method = "constant", ties = "ordered")
  return(list(logS=logS,logN=logN))
}

#
"piecewise.linear" <- function(x, y, nbreaks=1, middle=0.9){
  "piecewise.linear.likelihood" <- function(alpha, x, y) {
    nbreaks <- length(alpha)
    N <- length(x)
    w <- pmax(matrix(x,ncol=nbreaks,nrow=N) - matrix(alpha,ncol=nbreaks,nrow=N,byrow=T),0)
    X <- cbind(x,w)
    fit <- lm(y ~ X)
    SSE <- sum(fit$residuals^2)
    return(SSE)
  }
  r <- diff(range(x))
  offset <- r * (1 - middle)/2
  low <- min(x) + offset
  high <- max(x) - offset
  par0 <- as.vector(quantile(x,prob=seq((1-middle)/2,1-(1-middle)/2,len=nbreaks+2)[-c(1,nbreaks+2)]))
  temp <- optim(par0,fn=piecewise.linear.likelihood,x=x,y=y,method="L-BFGS-B",lower=low,upper=high)
  # final fit
  alpha <- temp$par
  N <- length(x)
  w <- pmax(matrix(x,ncol=nbreaks,nrow=N) - matrix(alpha,ncol=nbreaks,nrow=N,byrow=T),0)
  X <- cbind(x,w)
  fit <- lm(y ~ X)
  coeffs <- as.vector(coef(fit))
  #log10taus <- c(log10(exp(-coeffs[1]/coeffs[2])),alpha)
  log10taus <- c(low,alpha)
  betas <- cumsum(coeffs[-1])

  return(list(log10taus=log10taus,betas=betas))
}

#
initial <- function(Y, b, A, nbreaks){
  mest <- (Y-b)/A
  min0 <- min(mest[mest>0])
  mest[mest<=0] <- min0*1e-5
  temp <- logNlogS(mest)
  if (nbreaks==0){
    middle <- 0.9
    temp2 <- lm(temp$logN~temp$logS)
    betas <- -as.vector(coef(temp2))[2]
    r <- diff(range(temp$logS))
    offset <- r * (1 - middle)/2
    taus <- 10^(min(temp$logS) + offset)
  } else {
    temp2 <- piecewise.linear(x=temp$logS,y=temp$logN,nbreaks=nbreaks)
    taus.order <- order(temp2$log10taus)
    taus <- 10^temp2$log10taus[taus.order]
    betas <- -temp2$betas[taus.order]
  }
  return(list(taus=taus,betas=betas))
}

###########################
#### internal function ####
###########################

"Q3" <- function(tausm1,sam,const,useC=TRUE){
    Nsim <- nrow(sam)
    S <- array(sam)
    taus <- c(min(S),tausm1/const)
    taus <- sort(taus)
    optim.betas.res <- optim.betas(taus=taus,S=S,useC=useC)
    temp <- sum(log(optim.betas.res$betas)*optim.betas.res$mk)/Nsim
    if (is.na(temp)||is.infinite(temp)){
      return(9e300)
    }else{
      return(-temp)
    }
}

#########################################################
#### fit with chosen number of piece and (nonzero) b ####
#########################################################

"broken.power.B" <- function(Y, A, B,taus0=NULL,betas0=NULL,b=rep(0,length(Y)),Nlim=30,Nsim=7000,Nburn=2000,thin=5,plot.check=F,display=T,NM.maxit=200,NM.funevals=300,useC=TRUE){

  if (is.null(taus0)|is.null(betas0)){
    initial.res <- initial(Y,b,A,B-1) #note nbreaks = B-1
    taus0 <- initial.res$taus
    betas0 <- initial.res$betas
  }

  # checking
  if ((length(taus0)!=B)||(length(betas0)!=B)){
    stop("taus0 and betas0 are not consistent with B!\n")
  }

  # initialization
  Ncount <- 1
  diff.par <- 1
  pars.out1 <- matrix(nrow=Nlim+1,ncol=B*2)
  pars.out1[1,] <- c(taus0,betas0)
  pars.out2 <- pars.out1
  taus <- taus0
  betas <- betas0
  const <- mean(A)
  ES <- pmax((Y-b)/A,min(taus0))
  EU <- pBpareto(ES,taus0,betas0)

  # initialization for NelderMead optimization
  "Q3.NM" <- function(x=NULL, index=NULL, fmsfundata=NULL, useC=TRUE){
  	if (useC){
    	f <- Q3(as.vector(x),fmsfundata$sam,fmsfundata$const, useC=TRUE)
    } else {
		f <- Q3(as.vector(x),fmsfundata$sam,fmsfundata$const, useC=FALSE)
    }
    return(list(f=f, g=c(), c=c(), gc=c(), index=index, this=list(costfargument=fmsfundata)))
  }
  "Q.alt.NM" <- function(x=NULL, index=NULL, fmsfundata=NULL, useC=TRUE){
  	if (useC){
  		f <- Q.alt.b(as.vector(x),fmsfundata$Y,fmsfundata$A,fmsfundata$b,fmsfundata$sam,fmsfundata$const,useC=TRUE)
  	} else {
  		f <- Q.alt.b(as.vector(x),fmsfundata$Y,fmsfundata$A,fmsfundata$b,fmsfundata$sam,fmsfundata$const,useC=FALSE)
  	}
    return(list(f=f, g=c(), c=c(), gc=c(), index=index, this=list(costfargument=fmsfundata)))
  }

  nm1 <- neldermead.new()
  nm1 <- neldermead.configure(nm1,"-numberofvariables",B-1)
  nm1 <- neldermead.configure(nm1,"-function",Q3.NM)
  nm1 <- neldermead.configure(nm1,"-verbose",0)
  nm1 <- neldermead.configure(nm1,"-storehistory",F)
  nm1 <- neldermead.configure(nm1,"-verbosetermination",0)
  nm1 <- neldermead.configure(nm1,"-method","box")
  nm1 <- neldermead.configure(nm1,"-maxiter",NM.maxit)
  nm1 <- neldermead.configure(nm1,"-maxfunevals",NM.funevals)

  nm2 <- neldermead.new()
  nm2 <- neldermead.configure(nm2,"-numberofvariables",2*B)
  nm2 <- neldermead.configure(nm2,"-function",Q.alt.NM)
  nm2 <- neldermead.configure(nm2,"-verbose",0)
  nm2 <- neldermead.configure(nm2,"-storehistory",F)
  nm2 <- neldermead.configure(nm2,"-verbosetermination",0)
  nm2 <- neldermead.configure(nm2,"-method","box")
  nm2 <- neldermead.configure(nm2,"-maxiter",NM.maxit)
  nm2 <- neldermead.configure(nm2,"-maxfunevals",NM.funevals)


  # EM algorithm
  while ((diff.par>1e-10)&&(Ncount<=Nlim)){
    if (display){
      cat("###### Step:",Ncount,"#####\n")
    }

    ########################################################
    # Sufficient Data Augmentation
    if (display){
      cat("Sufficient Data Augmentation Scheme...\n","   E Step...\n")
    }

    # save
    taus.old <- taus
    betas.old <- betas
    ES.old <- ES

    # E step
    sam <- sample.post.S.b(Y=Y,A=A,taus=taus.old,betas=betas.old,b=b,S.old=ES.old,Nsim=Nsim,display=display,useC=useC)
    sam <- mcmc(sam)
    sam <- window(sam,start=Nburn+1,thin=thin)

    ## plot check
    if (plot.check){
      check.select <- sample(1:length(Y),10,replace=F)
      par(mfrow=c(5,2),mar=c(3,2,3,1))
      for (i in (1:10)){
       plot(as.matrix(sam)[,check.select[i]],main=paste(c("The ",check.select[i],"-th element of S"),collapse=""),type="l",xlab="",ylab="")
      }
    }

    # M step
    if (display){
      cat("   M Step...\n")
    }

    if (B==1){
      taus <- min(sam)
    }else{
      lower <- rep(min(sam)*const,B-1)
      upper <- rep(max(sam)*const,B-1)

      # NM optimization
      nm1 <- neldermead.configure(nm1,"-x0",transpose(pmin(pmax(taus.old[-1]*const,lower),upper)))
      nm1 <- neldermead.configure(nm1,"-boundsmin",lower)
      nm1 <- neldermead.configure(nm1,"-boundsmax",upper)
      fmsfundata <- list(sam=sam,const=const)
      attr(fmsfundata, "type") <- "T_FARGS"
      nm1 <- neldermead.configure(nm1,"-costfargument",fmsfundata)
      nm1.res <- neldermead.search(nm1)
      taus <- sort(c(min(sam),as.vector(neldermead.get(nm1.res,"-xopt"))/const))
    }

    betas <- optim.betas(taus=taus,S=array(sam),useC=useC)$betas
    pars.out1[Ncount+1,] <- c(taus,betas)

    ########################################################
    # Ancillary Data Augmentation
    if (display){
      cat("Ancillary Data Augmentation Scheme...\n","   E Step...\n")
    }

    # save
    taus.old <- taus
    betas.old <- betas
    EU.old <- EU

    # E step
    U.sam <- matrix(pBpareto(as.vector(sam),taus,betas),nrow=nrow(sam))
    ES <- apply(sam,2,mean)

    ## plot check
    if (plot.check){
      check.select <- sample(1:length(Y),10,replace=F)
      par(mfrow=c(5,2),mar=c(3,2,3,1))
      for (i in (1:10)){
	plot(as.matrix(sam)[,check.select[i]],main=paste(c("The ",check.select[i],"-th element of U"),collapse=""),type="l",xlab="",ylab="")
      }
    }

    # M step
    if (display){
      cat("   M Step...\n")
    }

    lower <- c(rep(1e-2,B),rep(0.01,B))
    upper <- c(rep(1e4,B),rep(10,B))

    # NM optimization
    nm2 <- neldermead.configure(nm2,"-x0",transpose(pmin(pmax(c(taus.old*const,betas.old),lower),upper)))
    nm2 <- neldermead.configure(nm2,"-boundsmin",lower)
    nm2 <- neldermead.configure(nm2,"-boundsmax",upper)
    fmsfundata <- list(Y=Y,A=A,b=b,sam=U.sam,const=const)
    attr(fmsfundata, "type") <- "T_FARGS"
    nm2 <- neldermead.configure(nm2,"-costfargument",fmsfundata)
    nm2.res <- neldermead.search(nm2)
    tpars <- as.vector(neldermead.get(nm2.res,"-xopt"))
    taus <- tpars[1:B]/const
    betas <- (tpars[-(1:B)])[order(taus)]
    taus <- sort(taus)

    ########################################################
    pars.out2[Ncount+1,] <- c(taus,betas)
    ########################################################

    #
    diff.par <- mean(((pars.out2[Ncount,]-pars.out2[Ncount+1,])/pars.out2[Ncount,])^2)
    if (display){
      cat("The difference is",diff.par,"and best parameter is",c(taus,betas),"\n\n")
    }
    Ncount <- Ncount+1
  }
  EU <- apply(sam,2,mean)
  return(list(pars=pars.out2[1:Ncount,],pars1=pars.out1[1:Ncount,],sam.U=sam,hat.pars=pars.out2[Ncount,], hat.taus=pars.out2[Ncount,1:B],hat.betas=pars.out2[Ncount,(B+1):(2*B)],EU=EU,ES=ES))
}

######################################################################################
#### compute loglikelihood (the log marginal likelihood) from the power posterior ####
######################################################################################

"loglike.power.post.b" <- function(Y,A,taus,betas,b,S.old,Nsim,Nburn,thin,c=3,Nint=50,display=T,multi=F,cpus=4,useC=TRUE){
  # the posterior expectation
  ## temperture schedule
  tmp <- seq(0,1,len=Nint)^c
  ## sampling form power posterior
  if (multi){
    registerDoMC(cores=cpus)
    temp <- foreach(mci = 1:Nint,.combine=c) %dopar% {
    	#function(Y,A,taus,betas,b,S.old,pwr,Nsim,display=T,useC=TRUE){
      sam <- window(mcmc(sample.power.post.S.b(Y=Y,A=A,taus=taus,betas=betas,b=b,S.old=S.old,pwr=tmp[mci],Nsim=Nsim,display=F,useC=useC)),start=Nburn+1,thin=thin)
      mean(apply(sam,1,function(x,Y,A,b){
      		sum(dpois(Y,lambda=(A*x+b),log=T))
      	},Y=Y,A=A,b=b))
    }
  } else {
    temp <- array(dim=Nint)
    for (k in (1:Nint)){
    	if (display){
    		cat("Sample from power posterior with temperture ", tmp[k],"...\n")
      	}
    	sam <- sample.power.post.S.b(Y=Y,A=A,taus=taus,betas=betas,b=b,S.old=S.old,pwr=tmp[k],Nsim=Nsim,display=F,useC=useC)
    	sam <- mcmc(sam)
    	sam <- window(sam,start=Nburn+1,thin=thin)
    	temp[k] <- mean(apply(sam,1,function(x,Y,A,b){sum(dpois(Y,lambda=(A*x+b),log=T))},Y=Y,A=A,b=b))
    	S.old <- apply(sam,2,mean)
    }
  }
  ## trapezoidal rule
  loglike <- trapezoidal(tmp,temp)
  return(loglike)
}

##########################################
#### model selection with (nonzero) b ####
##########################################
# fit a broken power law to the astronomical sources
"broken.power" <- function(Y,A,maxB=3,b=rep(0,length(Y)),Nlim=30,Nsim1=7000,Nburn1=2000,Nsim2=10000,Nburn2=2000,thin1=5,thin2=5,plot.check=F,display=F,PPmulti=T,cpus=4,NM.maxit=200,NM.funevals=300,details=F,useC=TRUE){
  est.pars <- array(list(),dim=maxB)
  IC <- matrix(nrow=maxB,ncol=3)
  bestB <- array(dim=3)
  loglike <- array(dim=maxB)
  out <- list(list())
  initial.res <- initial(Y, b, A, 0)
  taus0 <- initial.res$taus
  betas0 <- initial.res$betas

  for (j in (1:maxB)){
    cat("##################### B =",j," #####################\n")
    res <- broken.power.B(Y=Y,A=A,B=j,taus0=taus0,betas0=betas0,b=b,Nlim=Nlim,Nsim=Nsim1,Nburn=Nburn1,thin=thin1,plot.check=plot.check,display=display,NM.maxit=NM.maxit,NM.funevals=NM.funevals,useC=useC)
    est.pars[[j]] <- res$hat.pars
    if (display){
    	cat(paste0("Finished broken.power.B for B=",j,"...\n"))
    	cat("Calling loglike.power.post.b...\n")
    }
    #
    loglike[j] <- loglike.power.post.b(Y,A,res$hat.taus,res$hat.betas,b=b,res$ES,Nsim2,Nburn2,thin2,c=3,Nint=50,display=display,multi=PPmulti,cpus=cpus,useC=useC)
    if (display){
    	cat("done. Processing results...\n")
    }
    IC[j,1] <- -2*loglike[j] + 4*j
    IC[j,2] <- -2*loglike[j] + 2*log(length(Y))*j
    IC[j,3] <- -2*loglike[j] + 4*log(log(length(Y)))*j
    if (display){
    	cat("The estimated pars is ",est.pars[[j]],"\n")
    }
    initial.res <- initial(Y,b,A,j) #note nbreaks = B-1
    taus0 <- initial.res$taus
    betas0 <- initial.res$betas
    if (details) {
      out[[j]] <- res
    }
  }
  bestB <- apply(IC,2,which.min)
  B <- bestB[2]
  if (display){
  	cat("Best B's by criteria:\n") ; print(bestB)
  	cat("Best value for B:\n") ; print(B)
  	cat("Est.pars:\n") ; print(est.pars)
  }
  return(list(B=B,hat.taus=est.pars[[B]][1:B],hat.betas=est.pars[[B]][(B+1):(2*B)], est.pars=est.pars,IC=IC,bestB=bestB,loglike=loglike,out=out))
  #return(list(B=B,hat.taus=est.pars[[B]][,1:B],hat.betas=est.pars[[B]][,(B+1):(2*B)], est.pars=est.pars,IC=IC,bestB=bestB,loglike=loglike,out=out))
}
