## Plot the results of the power analysis.
library(ggplot2)
library(RColorBrewer)

scale1 <- brewer.pal(3, "Set1")
    
seed <- "12345"

this.theme <- theme_bw()+theme(legend.position="none", panel.border=element_rect(colour=NA), panel.grid.major=element_line(size=1), panel.grid.minor=element_line(size=0))

power <- read.table(paste0("~/selection/analysis/v8/power/read_power2_scale_1_seed_",seed,"_minf_0.1_N_1000.txt"), as.is=TRUE, header=TRUE)
names(power) <- c("Tag", "Pops", "s", "Generations", "Power", "Eff.N.chr")
power$Generations <- as.factor(power$Generations)
power$Pops <- as.factor(power$Pops)
plt <- ggplot(power, aes(x=s, y=Power, col=Generations, lty=Pops))+geom_line(lwd=1.5)+geom_point(pch=1, col="black")
## plt <- plt+scale_x_log10(breaks=c(seq(0.001, 0.01, 0.001), seq(0.02, 0.1, 0.01)))
plt <- plt+scale_x_log10(breaks=c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1))
plt <- plt+this.theme
pdf("~/selection/analysis/v8/power/read_power2.pdf")
print(plt)
dev.off()

power <- power[power$Pops=="All",]
plt <- ggplot(power, aes(x=s, y=Power, col=Generations))+geom_line(lwd=1.5)+geom_point(pch=1, col="black")
plt <- plt+scale_x_log10(breaks=c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1))
plt <- plt+this.theme+scale_color_manual(values=scale1)
pdf("~/selection/analysis/v8/power/read_power_all2.pdf")
print(plt)
dev.off()

## Now plot scale plot

group <- 2
scales <- c(1, 2, 5, 10)
ss <- c(0.005, 0.01, 0.02)
g <- 100
wh <- "All"
s <- 0.02

scale2 <- brewer.pal(length(scales), "Set2")
scale2 <- c(scale2[1], scale1[2], scale2[2:length(scale2)])

data <- read.table(paste0("~/selection/analysis/v8/power/read_power1_scale_1_seed_12345_minf_0.1_N_1000.txt"), as.is=TRUE, header=TRUE)
data <- data[data[,2]==wh & data[,4]==g,]
data$size=as.character(0.5)
results <- data

for(i in 1:length(scales)){
    data <- read.table(paste0("~/selection/analysis/v8/power/read_power2_scale_", scales[i], "_seed_12345_minf_0.1_N_1000.txt"), as.is=TRUE, header=TRUE)
    data <- data[data[,2]==wh & data[,4]==g,]
    data$size=as.character(scales[i])
    results <- rbind(results, data)
}

names(results) <- c("Tag", "pops", "s", "Generations", "Power", "Eff.N", "size")

this.theme <- theme_bw()+theme(legend.position="none", panel.border=element_rect(colour=NA), panel.grid.major=element_line(size=1), panel.grid.minor=element_line(size=0))
plt <- ggplot(results, aes(x=s, y=Power, col=size))+geom_line(lwd=1.5)+geom_point(pch=1, col="black")
plt <- plt+scale_x_log10(breaks=c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1))
plt <- plt+this.theme+scale_color_manual(values=scale2)


pdf("~/selection/analysis/v8/power/read_power_scale_all2.pdf")
print(plt)
dev.off()
