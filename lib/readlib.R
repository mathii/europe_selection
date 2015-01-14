## Functions for reading in data - mostly taking input from spindrift.

################################################################################
## Read in snp freq/count data generated with spinfridt/Freq.py
## root: file root
################################################################################

read.freq <- function(root){
    freq<-read.table(paste0(root,".freq"), as.is=TRUE, header=TRUE, na.strings="nan")
    count<-read.table(paste0(root,".count"), as.is=TRUE, header=TRUE, na.strings="nan")
    total<-read.table(paste0(root,".total"), as.is=TRUE, header=TRUE, na.strings="nan")
    freq <- convert(freq)
    count <- convert(count)
    total <- convert(total)
    return(list("freq"=freq, "count"=count, "total"=total))
}

################################################################################
## Combert spindrift data to a more useful format. 
################################################################################

convert <- function(freq){
    freq <- t(freq)
    colnames(freq) <- freq[1,]
    freq <- freq[6:NROW(freq),]
    mode(freq) <- "numeric"
    return(freq)
}


################################################################################
## Read in snp count/total data generated with spinfridt/Freq.py
## root: file root
################################################################################

read.counts.and.data <- function(root){    
    counts <- read.table(paste0(root, ".count"), header=TRUE, as.is=TRUE)
    totals <- read.table(paste0(root, ".total"), header=TRUE, as.is=TRUE)
    data <- counts[,1:5]
    counts <- data.matrix(counts[,6:NCOL(counts)])
    totals <- data.matrix(totals[,6:NCOL(totals)])
    return(list(data=data, counts=counts, totals=totals))
}
