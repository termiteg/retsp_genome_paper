d <- read.delim("R.B.biased.up.fdr0.01.genelist.gomat.BP.txt", comment.char="#", head=F, row.names=1)
colnames(d) <- c("target", "bg", "target.total", "bg.total", "target.percent", "bg.percent")
pval <- apply(d[,1:2], 1, 
  function(x){(fisher.test(matrix(c(x[1], 274-x[1], x[2], 13299-x[2]), 2), alternative="greater" ) )$p.value})
qval <- p.adjust(pval)
d <- cbind(d, pvalue=pval, FDR=qval)
write.table(d, "R.B.biased.up.fdr0.01.genelist.gomat.BP.fisher.Rout", sep="	", quote=F)
