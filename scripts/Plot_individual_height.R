## Plot the individual height loci for the German data.
## ~/spindrift/QPheno.py -d ~/data/v6/use/v61kg_europe2names -o indiv_height -p ~/data/v6/use/used_ancient_pseudodiploid.txt -i ~/data/v6/use/used_ancient_pseudodiploid.txt -e ~/selection/data/gwas/lango_allen_snps.gwas

morph.height <- read.table("~/selection/analysis/poly/Height/Nicole_morph_height", header=TRUE, as.is=TRUE, sep="\t")
morph.height$Individual <- gsub("^S", "I", morph.height$Individual)
genetic.height <- read.table("~/selection/analysis/poly/Height/QPheno_from_pulldown", header=TRUE, as.is=TRUE) 

data <- merge(genetic.height, morph.height, by="Individual")

data <- data[,c("Individual", "Population.x", "Infant", "Sex", "Morph.height", "Value")]
names(data)[2] <- "Population"
names(data)[6] <- "Genetic.height"
adult.data <- data[!data$Infant,]


pdf("~/selection/analysis/poly/Height/Qpheno_All.pdf")
pch.map=21:24
names(pch.map) <- c("LBK_EN", "Corded_Ware_LN", "Bell_Beaker_LN", "Unetice_EBA", "Alberstedt_LN")
coll <- ifelse(data$Infant, "darkred", "darkblue")
bgg <- ifelse(data$Infant, ifelse(data$Sex=="M", "darkred", "white"), ifelse(data$Sex=="M", "darkblue", "white") )
plot(data$Genetic.height, data$Morph.height, pch=pch.map[data$Population], col=coll, bg=bgg, bty="n",xlab="Genetic height", ylab="Morphological height")
legend("bottomleft", c(names(pch.map), "Male", "Female", "Infant"), pch=c(pch.map, 25, 25, 25), col=c(rep("darkblue", 6), "darkred"), pt.bg=c(rep("lightblue", 4), "darkblue", "white", "pink"), bty="n")
dev.off()

pdf("~/selection/analysis/poly/Height/Qpheno_Adult.pdf")
pch.map=21:25
names(pch.map) <- c("LBK_EN", "Corded_Ware_LN", "Bell_Beaker_LN", "Unetice_EBA", "Alberstedt_LN")
plot(adult.data$Genetic.height, adult.data$Morph.height, pch=pch.map[adult.data$Population], col="darkblue", bg=ifelse(adult.data$Sex=="M", "darkblue", "white"), bty="n", ylim=c(150,170), xlab="Genetic height", ylab="Morphological height", cex=2)
legend("bottomleft", c(names(pch.map), "Male", "Female"), pch=c(pch.map, 21, 21), col="darkblue", pt.bg=c(rep("lightblue", 5), "darkblue", "white"), bty="n")
dev.off()

