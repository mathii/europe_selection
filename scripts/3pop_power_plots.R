## Plot the results of the power analysis.
library(ggplot2)

this.theme <- theme_bw()+theme(legend.position="none", panel.border=element_rect(colour=NA), panel.grid.major=element_line(size=1), panel.grid.minor=element_line(size=0))

power <- read.table("~/selection/analysis/v8/power/read_power2.txt", as.is=TRUE, header=TRUE)
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
## plt <- plt+scale_x_log10(breaks=c(seq(0.001, 0.01, 0.001), seq(0.02, 0.1, 0.01)))
plt <- plt+scale_x_log10(breaks=c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1))
plt <- plt+this.theme
pdf("~/selection/analysis/v8/power/read_power_all2.pdf")
print(plt)
dev.off()
