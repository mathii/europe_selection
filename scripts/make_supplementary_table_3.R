## Read and merge the selection scan results with the frequency estimates.

f.order <- c("HG", "EF", "SA", "CEU", "GBR", "IBS", "TSI")

gscan <- read.table("~/selection/analysis/v8/gscan/scan_results_read2.corrected.txt", as.is=TRUE, header=TRUE)
freq <- read.table("~/selection/counts/v8/all.reads.3pop.freq", as.is=TRUE, header=TRUE)

mrg <- merge(gscan, freq, by="ID")

mrg <- mrg[,c("ID", "CHR", "POS", "REF", "ALT", "corrected.p", f.order)]
mrg[,f.order] <- round(1-mrg[,f.order],3)
mrg <- mrg[order(mrg$CHR, mrg$POS),]

mrg$corrected.p <- format(mrg$corrected.p, digits=3)

write.table(mrg, "~/selection/paper/spreadsheets/Supplementary_table_3.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
