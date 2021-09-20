library(bsseq)

d <- read.bismark("bismark_GenomeWide_Cytosine_Report_CpG_merged2.txt", "total", fileType="cytosineReport")
d.sm <- BSmooth(d, verbose=TRUE)

## data.frame to dump data
x <- as.data.frame(granges(d.sm))
x <- cbind(x, MethRaw=getMeth(d.sm, type="raw"))
x <- cbind(x, MethSmooth=getMeth(d.sm, type="smooth"))
colnames(x)[6:7] <- c("methylation_ratio_raw", "methylation_ratio_smoothed")
write.table(x, "methylation_ratio_bsmooth.txt", sep="\t", quote=F, row.names=F)