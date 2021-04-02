library(edgeR)

d <- read.delim("../fcounts/counts.txt", header=T, row.names=1)
colnames(d) <- toupper( sub("...tophat.", "", colnames(d)) )

len <- read.delim("../RsGM8.len.txt", header=F, row.names=1)
rownames(len) <- sub("-PA", "", rownames(len))


### Calc cpm & rpkm of all genes ###

y <- DGEList(counts=d)
y <- calcNormFactors(y, method="TMM")

cpm <- cpm(y)
len.ordered <- len[rownames(y$counts),]
rpkm <- rpkm(y, len.ordered)

write.table(cpm, file="cpm.tmm.allGenes.txt", sep="\t", quote=F)
write.table(rpkm, file="rpkm.tmm.allGenes.txt", sep="\t", quote=F)



### GLM analysis ###

Sample <- colnames(d)
Caste <- c(rep("Royal",12), rep("Soldier",12), rep("Worker",12))
Sex <- rep(c(rep("Female",6), rep("Male",6)), 3)
Part <- rep(c(rep("Body",3), rep("Head",3)), 6)

targets <- data.frame(Sample, Caste, Sex, Part)
design <- model.matrix(~ Part + Caste:Part + Sex:Part, data=targets)

y <- DGEList(counts=d)
keep <- rowSums(cpm(y)>=1) >= 3
y <- y[keep,, keep.lib.sizes=FALSE]

y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design, verbose=T)
# Disp = 0.04302 , BCV = 0.2074
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
colnames(fit)
#  [1] "(Intercept)"                   "PartHead"
#  [3] "PartBody:CasteSoldier"         "PartHead:CasteSoldier"
#  [5] "PartBody:CasteWorker"          "PartHead:CasteWorker"
#  [7] "PartBody:SexMale"              "PartHead:SexMale"


lrt_c.h <- glmLRT(fit, coef=c(4,6))
lrt_c.b <- glmLRT(fit, coef=c(3,5))
lrt_s.h <- glmLRT(fit, coef=8)
lrt_s.b <- glmLRT(fit, coef=7)
lrt_rs.h <- glmLRT(fit, coef=4)
lrt_rs.b <- glmLRT(fit, coef=3)
lrt_rw.h <- glmLRT(fit, coef=6)
lrt_rw.b <- glmLRT(fit, coef=5)

tt_c.h <- topTags(lrt_c.h, n=nrow(y$counts), adjust.method="BH")
tt_c.b <- topTags(lrt_c.b, n=nrow(y$counts), adjust.method="BH")
tt_s.h <- topTags(lrt_s.h, n=nrow(y$counts), adjust.method="BH")
tt_s.b <- topTags(lrt_s.b, n=nrow(y$counts), adjust.method="BH")
tt_rs.h <- topTags(lrt_rs.h, n=nrow(y$counts), adjust.method="BH")
tt_rs.b <- topTags(lrt_rs.b, n=nrow(y$counts), adjust.method="BH")
tt_rw.h <- topTags(lrt_rw.h, n=nrow(y$counts), adjust.method="BH")
tt_rw.b <- topTags(lrt_rw.b, n=nrow(y$counts), adjust.method="BH")

write.table(tt_c.h$table, file="LRT_Caste_Head.txt", sep="\t", quote=F)
write.table(tt_c.b$table, file="LRT_Caste_ThAb.txt", sep="\t", quote=F)
write.table(tt_s.h$table, file="LRT_Sex_Head.txt", sep="\t", quote=F)
write.table(tt_s.b$table, file="LRT_Sex_ThAb.txt", sep="\t", quote=F)
write.table(tt_rs.h$table, file="LRT_RvsS_Head.txt", sep="\t", quote=F)
write.table(tt_rs.b$table, file="LRT_RvsS_ThAb.txt", sep="\t", quote=F)
write.table(tt_rw.h$table, file="LRT_RvsW_Head.txt", sep="\t", quote=F)
write.table(tt_rw.b$table, file="LRT_RvsW_ThAb.txt", sep="\t", quote=F)

targets$Caste <- relevel(targets$Caste, ref="Soldier")
design <- model.matrix(~ Part + Caste:Part + Sex:Part, data=targets)
fit <- glmFit(y, design)
colnames(fit)
# [1] "(Intercept)"          "PartHead"
# [3] "PartBody:CasteRoyal"  "PartHead:CasteRoyal"
# [5] "PartBody:CasteWorker" "PartHead:CasteWorker"
# [7] "PartBody:SexMale"     "PartHead:SexMale"
lrt_sw.h <- glmLRT(fit, coef=6)
lrt_sw.b <- glmLRT(fit, coef=5)
tt_sw.h <- topTags(lrt_sw.h, n=nrow(y$counts), adjust.method="BH")
tt_sw.b <- topTags(lrt_sw.b, n=nrow(y$counts), adjust.method="BH")
write.table(tt_sw.h$table, file="LRT_SvsW_Head.txt", sep="\t", quote=F)
write.table(tt_sw.b$table, file="LRT_SvsW_ThAb.txt", sep="\t", quote=F)



### calc cpm & rpkm of statistically tested genes ###
 
cpm <- cpm(y)
write.table(cpm, file="cpm.tmm.glm.txt", sep="\t", quote=F)
 
len.keep <- len[rownames(y$counts),]
rpkm <- rpkm(y, len.keep)
write.table(cpm, file="rpkm.tmm.glm.txt", sep="\t", quote=F)