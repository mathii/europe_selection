#Functions to implement the 3 (or n) population test.
## require(numDeriv)

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

#########################################################
#
# Likelihood for the combined read and count based
# model.
# p: Reference allele frequencies (length N)
# data: list of N populations, each entry is a list with
# two elements, "reads" and "counts". "reads" is itself
# a list with "ref" and "alt" vectors of ref and alt
# read counts for each individual and "counts" is a vector
# of length 2 for giving ref/alt counts.
# This returns the log likelihood summed over all the
# populations. 
#########################################################

likelihood.reads <- function(freq, data, error.prob=0.001){
    ll <- 0
    i.pop <- 1
    for(pop in names(data)){
        N.read.ind <- length(data[[pop]][["reads"]][["ref"]]) #Number of individuals with read information
        p <- freq[i.pop]
        if(N.read.ind){
            #Add reads
            for(i in 1:N.read.ind){
                ref <- data[[pop]][["reads"]][["ref"]][[i]]
                alt <- data[[pop]][["reads"]][["alt"]][[i]]
                ## ll <- ll+log((alt==0)*p*p + (ref==0)*(1-p)*(1-p) + 2*dbinom(ref, ref+alt, 0.5)*p*(1-p))
                ll <- ll+log(dbinom(alt, ref+alt, error.prob)*p*p + dbinom(ref, ref+alt, error.prob)*(1-p)*(1-p) + dbinom(ref, ref+alt, 0.5)*2*p*(1-p))
            }
        }
        #Add counts, can be zero
        if("counts" %in% names(data[[pop]])){
            cnts <- data[[pop]][["counts"]]
            ll <- ll+dbinom(cnts[1], cnts[1]+cnts[2], p, log=TRUE)
        }
        i.pop <- i.pop+1    
    }
    return(ll)
}

#########################################################
#
# Constrained likelihood for read based model
# 
#########################################################

constrained.likelihood.reads <- function(p, A, data, error.prob=0){
    p <- c(p, c(p %*% A))
    p <- pmin(pmax(EPSILON, p), 1-EPSILON)
    ll <- likelihood.reads(p, data, error.prob=error.prob)
    return(ll)
}

#########################################################
#
# Fit the constrained model
# 
#########################################################

fit.constrained.model.reads <- function(data, A, error.prob=0){
    p.anc.init <- rep(0.5, dim(A)[1])
    opt <- optim(p.anc.init, constrained.likelihood.reads, A=A, data=data, error.prob=error.prob, control=list(fnscale=-1), lower=0, upper=1, method="L-BFGS")
    return(opt)
}

#########################################################
#
# Fit the unconstrained model
# Since the frequencies are independent, we can fit them
# separately. 
# 
#########################################################

fit.unconstrained.model.reads <- function(data, error.prob=0){
    p <- rep(NA, length(data))
    for(i in 1:length(data)){
        subdata=list(data[[i]])
        names(subdata) <- names(data)[i]
        opt <- optimize(likelihood.reads, data=subdata, error.prob=error.prob, maximum=TRUE, lower=EPSILON, upper=1-EPSILON)
        p[i] <- opt$maximum
    }
    ll <- likelihood.reads(p, data)
    return(list(par=p, value=ll))
}

#########################################################
#
# General 3 population test using read information                                        
# 
#########################################################

test.3pop.reads <- function(data, A, error.prob=error.prob){
    degf <- dim(A)[2]

    if(sum(dim(A))!=length(data)){stop("Matrix A not compatible with observations")}

    uf <- fit.unconstrained.model.reads(data, error.prob=error.prob)
    cf <- fit.constrained.model.reads(data, A, error.prob=error.prob)   
    
    stat <- 2*(uf$value-cf$value)
    p <- pchisq(stat, df=degf, lower.tail=F)
    return( c(stat, p) )
}

#########################################################
#
# Estimated frequency and effective sample size for
# one data point. Freq.data of size 1 
#########################################################

effective.size.reads <- function(data, error.prob=0){
    uf <- fit.unconstrained.model.reads(data, error.prob=error.prob)
    d2 <- hessian(likelihood.reads, uf$par, method="Richardson", method.args=list(eps=EPSILON/2), data, error.prob=error.prob)
    N <- -uf$par*(1-uf$par)*c(d2)
    return(list(p=uf$par, N=N, d2=c(d2)))
}

#########################################################
#
# Empty read/allele count data structure
#
#########################################################

make.empty.data <- function(pops){
    empty.data <- rep( list(list()), length(pops) ) 
    names(empty.data) <- pops
    for(pop in pops){
        empty.data[[pop]] <- list("reads"=list("ref"=NULL, "alt"=NULL), "counts"=c(0,0)) #ref and alt counts. 
    }
    return(empty.data)
}

#########################################################
#
# Combine information from read counts and allele counts 
# into the data structure used for the likelihood
# calcluations here
#
# pops: list of populations to include
# include.reads: list, with names in pops, mapping to
# vectors of sub-populations to add together to make
# each population, using read counts.
# include.read samples: which samples are in each sub-population
# include.counds: As include.reads, but with populations
# for which we should add hard counts.
# read.info, counts, totals: standard data structures
#########################################################

make.freq.data <- function(pops, include.reads, include.read.samples, include.counts, read.info, counts, totals, empty.data){
    freq.data <- empty.data
    for(pop in pops){
        if(pop %in% names(include.reads)){
            for(sample in include.read.samples[[pop]]){
                ref.alt <- read.info[read.info[,2]==sample,3:4]
                if(sum(ref.alt)>0){
                    freq.data[[pop]][["reads"]][["ref"]] <- c(freq.data[[pop]][["reads"]][["ref"]],ref.alt[[1]])
                    freq.data[[pop]][["reads"]][["alt"]] <- c(freq.data[[pop]][["reads"]][["alt"]],ref.alt[[2]])
                }
            }
        }

        if(pop %in% names(include.counts)){
            for(subpop in include.counts[[pop]]){
                ref.alt <- c(counts[i,subpop], totals[i,subpop]-counts[i,subpop])
                names(ref.alt) <- NULL
                freq.data[[pop]][["counts"]] <- freq.data[[pop]][["counts"]]+ref.alt
            }
        }
    }
    return(freq.data)
}
