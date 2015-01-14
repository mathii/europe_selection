#Functions to implement the 3 (or n) population test.

EPSILON <- 1e-4 #Minumum allele frequency

#########################################################
#
# Likelihood of allele counts N.A/N given population
# frequencies p in source and Ap in derived pops
#
#########################################################

constrained.likelihood <- function(p, A, N, N.A){
    p <- c(p, c(p %*% A))
    p <- pmin(pmax(EPSILON, p), 1-EPSILON)
    return(sum(dbinom(N.A, N, p, log=TRUE)))
}

#########################################################
#
# Fit the allele frequencies in the source populations. 
#                                        
#########################################################


fit.constrained.model <- function(N, N.A, A){
    p.anc.init <- rep(0.5, dim(A)[1])
    opt <- optim(p.anc.init, constrained.likelihood, A=A, N=N, N.A=N.A, control=list(fnscale=-1), lower=0, upper=1, method="L-BFGS")
    return(opt)
}

#########################################################
#
# Likelihood ratio test, given observations N.A/N and
# matrix A mapping source to derived populations.
# N and N.A are length K and A is a matrix with
# sum(dim(A))==K
#
#########################################################

test.3pop <- function(N, N.A, A){
    degf <- dim(A)[2]

    if(length(N)!=length(N.A)){stop("N and N.A different lengths")}
    if(sum(dim(A))!=length(N)){stop("Matrix A not compatible with observations")}

    #Alternative
    #MLE is just the observed proportion
    #If no observations then it doesn't matter what p is
    p.hat <- ifelse(N==0, 0.5, N.A/N)  
    p.hat <- pmin(pmax(EPSILON, p.hat), 1-EPSILON)
    l1 <- sum(dbinom(N.A, N, p.hat, log=TRUE))

    #Null
    opt <- fit.constrained.model(N, N.A, A)
    p.est <- c(opt$par, opt$par %*% A) 
    p.est <- pmin(pmax(EPSILON, p.est), 1-EPSILON)
    l0 <- sum(dbinom(N.A, N, p.est, log=TRUE))

    stat <- 2*(l1-l0)
    p <- pchisq(stat, df=degf, lower.tail=F)
    return( c(stat, p) )
}

#########################################################
#
# Likelihood ratio test, given observations N.A/N, test
# whether all the observations can come from the
# same population. 
#
#########################################################

test.diff <- function(N, N.A){
    degf <- length(N)-1

    if(length(N)!=length(N.A)){stop("N and N.A different lengths")}

    #Alternative
    #MLE is just the observed proportion
    #If no observations then it doesn't matter what p is
    p.hat <- ifelse(N==0, 0.5, N.A/N)  
    p.hat <- pmin(pmax(EPSILON, p.hat), 1-EPSILON)
    l1 <- sum(dbinom(N.A, N, p.hat, log=TRUE))

    #Null
    p.est <- rep(sum(N.A)/sum(N), length(N))
    p.est <- pmin(pmax(EPSILON, p.est), 1-EPSILON)
    l0 <- sum(dbinom(N.A, N, p.est, log=TRUE))

    stat <- 2*(l1-l0)
    p <- pchisq(stat, df=degf, lower.tail=F)
    return( c(stat, p) )
}
