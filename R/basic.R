####################################
#### pdf of pareto distribution ####
####################################

dpareto <- function(x,tau,beta){
    temp <- beta*tau^beta/x^(beta+1)
    ifelse(x>=tau,temp,0)
}

######################################################
#### random generator of a pareto random variable ####
######################################################

rpareto <- function(n,tau,beta){
    tau*(runif(n))^(-1/beta)
}

#######################################
#### log incomplete gamma function ####
#######################################
# require gsl library

# logigamma <- function(a,x){
#     # int^\infty_x t^{a-1} e^{-t} dt
#     res <- array(dim=length(a))
#     ind <- (a<=0)
#     ind1 <- !ind
#     if (sum(ind)>0){
# 	#cat(a[ind],"\n")
# 	temp1 <- log(gamma_inc(a[ind],x[ind]))
#         res[ind] <- temp1
#     }
#     if (sum(ind1)>0){
# 	temp2 <- pgamma(x[ind1],a[ind1],lower.tail=F,log.p=T)+lgamma(a[ind1])
# 	res[ind1] <- temp2
#     }
#     return(res)
# }

########################################################################
#### negative log likelihood of a single piece poisson pareto model ####
########################################################################

# neg.log.like.single <- function(alpha,Y,A,m){
#     #cat(alpha,m,"\n")
#     # alpha: the exponent parameter of the pareto distribution
#     # m: the min. parameter of the pareto distribution
#     return(-(length(Y)*log(alpha) + alpha*sum(log(m*A)) - sum(lgamma(Y+1)) + sum(logigamma(Y-alpha,A*m))))
# }


#####################################################
#### MLE for a single piece poisson pareto model ####
#####################################################

# wrapper function for pois.pareto1

# neg.log.like.single.wrapper <- function(par,Y,A){
#     neg.log.like.single(par[2],Y,A,par[1]) #changed
# }

# pois.pareto1 <- function(Y,A,par0=c(1,1)){
#     const <- mean(A)
#     tempA <- A/const
#     result <- optim(par=par0,fn=neg.log.like.single.wrapper,Y=Y,A=tempA,method="L-BFGS-B",lower=c(1e-5/max(tempA),0.01),upper=c(20,20))
#     est <- c(result$par[1]/const,result$par[2])
#     return(list(est=est,nllike=result$value,optim.res=result))
# }

###############################
#### pdf of B piece pareto ####
###############################
#
dBpareto <- function(x,taus,betas){
    # taus and betas are of the same length
    betas <- betas[order(taus)]
    taus <- sort(taus)
    B <- length(taus)
    const <- c(1,cumprod((taus[-B]/taus[-1])^betas[-B]))
    temp <- (x>=taus[1])*1
    taus1 <- c(taus,Inf)
    for (j in (1:B)){
	temp <- ifelse(((x>=taus1[j])*(x<taus1[j+1]))==1,const[j]*betas[j]*taus[j]^betas[j]/x^(betas[j]+1),temp)
    }
    return(temp)
}

#
log.dBpareto <- function(x,taus,betas){
    # taus and betas are of the same length
    betas <- betas[order(taus)]
    taus <- sort(taus)
    B <- length(taus)
    const <- c(0,cumsum(betas[-B]*(log(taus[-B])-log(taus[-1]))))
    #const <- c(1,cumprod((taus[-B]/taus[-1])^betas[-B]))
    temp <- (x<taus[1])*(-Inf)
    taus1 <- c(taus,Inf)
    for (j in (1:B)){
	temp <- ifelse(((x>=taus1[j])*(x<taus1[j+1]))==1,const[j] + log(betas[j])+log(taus[j])*betas[j]-log(x)*(betas[j]+1),temp)
    }
    return(temp)
}

#
pBpareto <- function(x,taus,betas){
    B <- length(taus)
    betas <- betas[order(taus)]
    taus <- sort(taus)
    const <- c(1,cumprod((taus[-B]/taus[-1])^betas[-B]))
    temp <- (x>=taus[1])*1
    taus1 <- c(taus,Inf)
    for (j in (1:B)){
	temp <- ifelse(((x>=taus1[j])*(x<taus1[j+1]))==1,1-const[j]*(taus[j]/x)^betas[j],temp)
    }
    return(temp)
}

#########################
#### data simulation ####
#########################

sim.dat.pareto.pois.b <- function(n,taus,betas,A,b,useC=TRUE){
    B <- length(taus)
    S <- rBpareto(n,taus,betas,useC=useC)
    Y <- rpois(n,A*S+b)
    return(list(Y=Y,A=A,S=S,b=b))
}

##########################
#### trapezoidal rule ####
##########################

trapezoidal <- function(x,y){
    n <- length(y)
    return(sum(diff(x)*(y[-n]+y[-1])/2))
}


