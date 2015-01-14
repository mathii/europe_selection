#Process raw v6 population names.
#1. Replace pop names with Iosif's names.
#2. Rename everything that's WHG to WHG.
pop.labels <- c("ESN", "GWD", "LWK", "MSL", "YRI", "ACB", "ASW", "CLM", "MXL", "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU")


indfile <- read.table("~/data/v6/raw/v61kgx.ind", header=FALSE, as.is=TRUE)
rownames(indfile) <- indfile[,1]

europe2names <- read.table("~/data/v6/raw/europe2.ind", header=FALSE, as.is=TRUE)

for(i in 1:NROW(europe2names)){
    if( europe2names[i,1] %in% rownames(indfile) & !(indfile[europe2names[i,1],3] %in% pop.labels)){
        indfile[europe2names[i,1],3] <- europe2names[i,3] 
    }
}

write.table(indfile, "~/data/v6/use/v61kg_europe2names.ind", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")


## LBK.EN.samples <- c()
## WHG.samples <- c("Loschbour", "I0572", "KO1", "LaBrana1")

## for(w in WHG.samples){
##     indfile[w,3] <- "WHG"
## }

## write.table(indfile, "~/data/v6/use/v61kg_europe2names_WHGmerge.ind", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
