library(ggplot2)
data <- read.delim("rpkm.tmm.allGenes.txt", header=T)
 
part <- rep( c(rep("thorax and abdomen", 3), rep("head", 3)), 6)
sex <- rep( c(rep("female", 6), rep("male", 6)), 3)
caste <- c(rep("royal", 12), rep("soldier", 12), rep("worker", 12))
 
for (i in 1:nrow(data)) {
 
  d <- data[i,]
  df <- data.frame("RPKM"=as.numeric(d), part, sex, caste)
  min <- min(d)
 
  pngname <- paste("plot_rpkm_ymin0/", rownames(d), ".CasteSex.rpkm.png", sep="")
  png(filename=pngname, width=1000, height=1400, pointsize=30, bg="white")
  g <- ggplot(df,aes(x=caste, y=RPKM, colour=sex, fill=sex))
  g <- g + ggtitle(rownames(d)) + facet_wrap(~part, ncol=1, scales="free_y") + geom_hline(yintercept=0, alpha=0)
  g <- g + geom_point(position=position_jitterdodge(0.25, dodge.width=0.5), size=18, alpha=0.65, shape=20)
  g <- g + theme(
    title=element_text(size=rel(4.5)),
    axis.title.x=element_text(size=rel(1.3), face="bold"),
    axis.title.y=element_text(size=rel(1.2), face="bold"),
    axis.text.x = element_text(size=rel(4.8), colour="#444444"),
    axis.text.y=element_text(size=rel(4.8), colour="#444444"),
    strip.text=element_text(size=rel(4), face="bold"),
    legend.title=element_blank(),
    legend.text=element_text(size=rel(3.5)),
    legend.key=element_blank()
    )
  plot(g)
  dev.off()
 
}
