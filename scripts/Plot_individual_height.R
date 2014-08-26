## Plot the individual height loci for the German data.

morph.height <- read.table("~/selection_map/height/Individual/Nicole_morph_height", header=TRUE, as.is=TRUE, sep="\t")
genetic.height <- read.table("~/selection_map/height/Individual/QPheno_from_pulldown", header=TRUE, as.is=TRUE) 

data <- merge(genetic.height, morph.height, by="Individual")

data <- data[,c("Individual", "Population.x", "Infant", "Sex", "Morph.height", "Value")]
names(data)[2] <- "Population"
names(data)[6] <- "Genetic.height"
adult.data <- data[!data$Infant,]


pdf("~/selection_map/height/Individual/All.pdf")
pch.map=21:24
names(pch.map) <- c("German_7000BP_LBK", "German_4500BP_CordedWare", "German_4200BP_BellBeaker", "German_4000BP_Unetice")
coll <- ifelse(data$Infant, "darkred", "darkblue")
bgg <- ifelse(data$Infant, ifelse(data$Sex=="M", "darkred", "white"), ifelse(data$Sex=="M", "darkblue", "white") )
plot(data$Genetic.height, data$Morph.height, pch=pch.map[data$Population], col=coll, bg=bgg, bty="n",xlab="Genetic height", ylab="Morphological height")
legend("bottomleft", c(names(pch.map), "Male", "Female", "Infant"), pch=c(pch.map, 25, 25, 25), col=c(rep("darkblue", 6), "darkred"), pt.bg=c(rep("lightblue", 4), "darkblue", "white", "pink"), bty="n")
dev.off()

pdf("~/selection_map/height/Individual/Adult.pdf")
pch.map=21:24
names(pch.map) <- c("German_7000BP_LBK", "German_4500BP_CordedWare", "German_4200BP_BellBeaker", "German_4000BP_Unetice")
plot(adult.data$Genetic.height, adult.data$Morph.height, pch=pch.map[adult.data$Population], col="darkblue", bg=ifelse(adult.data$Sex=="M", "darkblue", "white"), bty="n", ylim=c(150,170), xlab="Genetic height", ylab="Morphological height")
legend("bottomleft", c(names(pch.map), "Male", "Female"), pch=c(pch.map, 25, 25), col="darkblue", pt.bg=c(rep("lightblue", 4), "darkblue", "white"), bty="n")
dev.off()

