#Rename the populations as we want to use them for the
#polygenic analysis.
#Don't forget to create symlinks to the snp and geno files. 

data <- read.table("~/data/v6/use/v61kg_europe2names.ind", as.is=TRUE)
polymap <- read.table("~/data/v6/use/polymap.txt", as.is=TRUE)

for(i in 1:NROW(polymap)){
    data[data[,3]==polymap[i,1],3] <- polymap[i,2]
}

write.table(data, "~/data/v6/use/v61kg_polynames.ind", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
