

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
# design <- model.matrix(~ Part + Caste:Part + Sex:Part + Caste:Sex:Part, data=targets)
 
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
# lrt_p.r <- glmLRT(fit, contrast=c(0,0,-1,1,0,0,0,0))
tt_sw.h <- topTags(lrt_sw.h, n=nrow(y$counts), adjust.method="BH")
tt_sw.b <- topTags(lrt_sw.b, n=nrow(y$counts), adjust.method="BH")
# tt_p.r <- topTags(lrt_p.r, n=nrow(y$counts), adjust.method="BH")
write.table(tt_sw.h$table, file="LRT_SvsW_Head.txt", sep="\t", quote=F)
write.table(tt_sw.b$table, file="LRT_SvsW_ThAb.txt", sep="\t", quote=F)
# write.table(tt_p.s$table, file="LRT_Part_Royal.txt", sep="\t", quote=F)
 
 
 
### calc cpm & rpkm of statistically tested genes ###
 
cpm <- cpm(y)
write.table(cpm, file="cpm.tmm.glm.txt", sep="\t", quote=F)
 
len.keep <- len[rownames(y$counts),]
rpkm <- rpkm(y, len.keep)
write.table(cpm, file="rpkm.tmm.glm.txt", sep="\t", quote=F)


GLM, Pairwise comparisons for differences between sexes

library(edgeR)
 
d <- read.delim("../fcounts/counts.txt", header=T, row.names=1)
colnames(d) <- toupper( sub("...tophat.", "", colnames(d)) )
 
len <- read.delim("../RsGM8.len.txt", header=F, row.names=1)
rownames(len) <- sub("-PA", "", rownames(len))
 
 
### Pairwirse comparison between females and males of each caste ###
 
Sample <- colnames(d)
Caste <- c(rep("Royal",12), rep("Soldier",12), rep("Worker",12))
Sex <- rep(c(rep("Female",6), rep("Male",6)), 3)
Part <- rep(c(rep("Body",3), rep("Head",3)), 6)
 
Group <- paste(Caste, Sex, Part, sep=".")
 
targets <- data.frame(Group)
design <- model.matrix(~ 0 + Group, data=targets)
 
y <- DGEList(counts=d)
keep <- rowSums(cpm(y)>=1) >= 3
y <- y[keep,, keep.lib.sizes=FALSE]
 
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design, verbose=T) # => Disp = 0.04201 , BCV = 0.205
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
 
fit <- glmFit(y, design)
colnames(fit)
#  [1] "GroupRoyal.Female.Body"   "GroupRoyal.Female.Head"
#  [3] "GroupRoyal.Male.Body"     "GroupRoyal.Male.Head"
#  [5] "GroupSoldier.Female.Body" "GroupSoldier.Female.Head"
#  [7] "GroupSoldier.Male.Body"   "GroupSoldier.Male.Head"
#  [9] "GroupWorker.Female.Body"  "GroupWorker.Female.Head"
# [11] "GroupWorker.Male.Body"    "GroupWorker.Male.Head"
 
lrt_RmRf_b <- glmLRT(fit, contrast=c(1,0,-1,0,0,0,0,0,0,0,0,0))
lrt_RmRf_h <- glmLRT(fit, contrast=c(0,1,0,-1,0,0,0,0,0,0,0,0))
lrt_SmSf_b <- glmLRT(fit, contrast=c(0,0,0,0,1,0,-1,0,0,0,0,0))
lrt_SmSf_h <- glmLRT(fit, contrast=c(0,0,0,0,0,1,0,-1,0,0,0,0))
lrt_WmWf_b <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,1,0,-1,0))
lrt_WmWf_h <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,0,0,1,0,-1))
 
tt_RmRf_b <- topTags(lrt_RmRf_b, n=nrow(y$counts), adjust.method="BH")
tt_RmRf_h <- topTags(lrt_RmRf_h, n=nrow(y$counts), adjust.method="BH")
tt_SmSf_b <- topTags(lrt_SmSf_b, n=nrow(y$counts), adjust.method="BH")
tt_SmSf_h <- topTags(lrt_SmSf_h, n=nrow(y$counts), adjust.method="BH")
tt_WmWf_b <- topTags(lrt_WmWf_b, n=nrow(y$counts), adjust.method="BH")
tt_WmWf_h <- topTags(lrt_WmWf_h, n=nrow(y$counts), adjust.method="BH")
 
write.table(tt_RmRf_b$table, file="LRT_RMvsRF_ThAb.txt", sep="\t", quote=F)
write.table(tt_RmRf_h$table, file="LRT_RMvsRF_Head.txt", sep="\t", quote=F)
write.table(tt_SmSf_b$table, file="LRT_SMvsSF_ThAb.txt", sep="\t", quote=F)
write.table(tt_SmSf_h$table, file="LRT_SMvsSF_Head.txt", sep="\t", quote=F)
write.table(tt_WmWf_b$table, file="LRT_WMvsWF_ThAb.txt", sep="\t", quote=F)
write.table(tt_WmWf_h$table, file="LRT_WMvsWF_Head.txt", sep="\t", quote=F)
 
 
 
### Pairwise comparison btwn a specific sex/caste and the others ###
 
my.contrasts.body <- makeContrasts(
  RfvsRm = GroupRoyal.Female.Body - GroupRoyal.Male.Body, 
  RfvsSf = GroupRoyal.Female.Body - GroupSoldier.Female.Body,
  RfvsSm = GroupRoyal.Female.Body - GroupSoldier.Male.Body,
  RfvsWf = GroupRoyal.Female.Body - GroupWorker.Female.Body,
  RfvsWm = GroupRoyal.Female.Body - GroupWorker.Male.Body,
  RmvsRf = GroupRoyal.Male.Body - GroupRoyal.Female.Body,
  RmvsSf = GroupRoyal.Male.Body - GroupSoldier.Female.Body,
  RmvsSm = GroupRoyal.Male.Body - GroupSoldier.Male.Body,
  RmvsWf = GroupRoyal.Male.Body - GroupWorker.Female.Body,
  RmvsWm = GroupRoyal.Male.Body - GroupWorker.Male.Body,
  SfvsRf = GroupSoldier.Female.Body - GroupRoyal.Female.Body,
  SfvsRm = GroupSoldier.Female.Body - GroupRoyal.Male.Body,
  SfvsSm = GroupSoldier.Female.Body - GroupSoldier.Male.Body,
  SfvsWf = GroupSoldier.Female.Body - GroupWorker.Female.Body,
  SfvsWm = GroupSoldier.Female.Body - GroupWorker.Male.Body,
  SmvsRf = GroupSoldier.Male.Body - GroupRoyal.Female.Body,
  SmvsRm = GroupSoldier.Male.Body - GroupRoyal.Male.Body,
  SmvsSf = GroupSoldier.Male.Body - GroupSoldier.Female.Body,
  SmvsWf = GroupSoldier.Male.Body - GroupWorker.Female.Body,
  SmvsWm = GroupSoldier.Male.Body - GroupWorker.Male.Body,
  WfvsRf = GroupWorker.Female.Body - GroupRoyal.Female.Body,
  WfvsRm = GroupWorker.Female.Body - GroupRoyal.Male.Body,
  WfvsSf = GroupWorker.Female.Body - GroupSoldier.Female.Body,
  WfvsSm = GroupWorker.Female.Body - GroupSoldier.Male.Body,
  WfvsWm = GroupWorker.Female.Body - GroupWorker.Male.Body,
  WmvsRf = GroupWorker.Male.Body - GroupRoyal.Female.Body,
  WmvsRm = GroupWorker.Male.Body - GroupRoyal.Male.Body,
  WmvsSf = GroupWorker.Male.Body - GroupSoldier.Female.Body,
  WmvsSm = GroupWorker.Male.Body - GroupSoldier.Male.Body,
  WmvsWf = GroupWorker.Male.Body - GroupWorker.Female.Body,
  levels = design
)
 
my.contrasts.head <- makeContrasts(
  RfvsRm = GroupRoyal.Female.Head - GroupRoyal.Male.Head, 
  RfvsSf = GroupRoyal.Female.Head - GroupSoldier.Female.Head,
  RfvsSm = GroupRoyal.Female.Head - GroupSoldier.Male.Head,
  RfvsWf = GroupRoyal.Female.Head - GroupWorker.Female.Head,
  RfvsWm = GroupRoyal.Female.Head - GroupWorker.Male.Head,
  RmvsRf = GroupRoyal.Male.Head - GroupRoyal.Female.Head,
  RmvsSf = GroupRoyal.Male.Head - GroupSoldier.Female.Head,
  RmvsSm = GroupRoyal.Male.Head - GroupSoldier.Male.Head,
  RmvsWf = GroupRoyal.Male.Head - GroupWorker.Female.Head,
  RmvsWm = GroupRoyal.Male.Head - GroupWorker.Male.Head,
  SfvsRf = GroupSoldier.Female.Head - GroupRoyal.Female.Head,
  SfvsRm = GroupSoldier.Female.Head - GroupRoyal.Male.Head,
  SfvsSm = GroupSoldier.Female.Head - GroupSoldier.Male.Head,
  SfvsWf = GroupSoldier.Female.Head - GroupWorker.Female.Head,
  SfvsWm = GroupSoldier.Female.Head - GroupWorker.Male.Head,
  SmvsRf = GroupSoldier.Male.Head - GroupRoyal.Female.Head,
  SmvsRm = GroupSoldier.Male.Head - GroupRoyal.Male.Head,
  SmvsSf = GroupSoldier.Male.Head - GroupSoldier.Female.Head,
  SmvsWf = GroupSoldier.Male.Head - GroupWorker.Female.Head,
  SmvsWm = GroupSoldier.Male.Head - GroupWorker.Male.Head,
  WfvsRf = GroupWorker.Female.Head - GroupRoyal.Female.Head,
  WfvsRm = GroupWorker.Female.Head - GroupRoyal.Male.Head,
  WfvsSf = GroupWorker.Female.Head - GroupSoldier.Female.Head,
  WfvsSm = GroupWorker.Female.Head - GroupSoldier.Male.Head,
  WfvsWm = GroupWorker.Female.Head - GroupWorker.Male.Head,
  WmvsRf = GroupWorker.Male.Head - GroupRoyal.Female.Head,
  WmvsRm = GroupWorker.Male.Head - GroupRoyal.Male.Head,
  WmvsSf = GroupWorker.Male.Head - GroupSoldier.Female.Head,
  WmvsSm = GroupWorker.Male.Head - GroupSoldier.Male.Head,
  WmvsWf = GroupWorker.Male.Head - GroupWorker.Female.Head,
  levels = design
)
 
 
spglist_body <- list() ; spglist_head <- list()
for ( i in 0:5 ) {
  specific_genes_b <- list() ; specific_genes_h <- list()
  for ( j in 1:5 ) {
    num <- j + i*5 ; print(num)
    lrt_b <- glmLRT(fit, contrast=my.contrasts.body[,num])
    lrt_h <- glmLRT(fit, contrast=my.contrasts.head[,num])
    tt_b <- topTags(lrt_b, n=nrow(y$counts), adjust.method="BH")
    tt_h <- topTags(lrt_h, n=nrow(y$counts), adjust.method="BH")
    sp_b <- rownames(subset(tt_b$table, logFC>3.32192809488736 & FDR<0.01))
    sp_h <- rownames(subset(tt_h$table, logFC>3.32192809488736 & FDR<0.01))
    specific_genes_b <- c(specific_genes_b, list(sp_b))
    specific_genes_h <- c(specific_genes_h, list(sp_h))
  }
  spglist_body <- c( spglist_body, list(Reduce(intersect, specific_genes_b)) )
  spglist_head <- c( spglist_head, list(Reduce(intersect, specific_genes_h)) )
}
 
listnames <- rownames(my.contrasts.body)[seq(1, nrow(my.contrasts.body), +2)]
listnames <- sub("Group", "", listnames)
listnames <- sub(".Body", "_specific_genes", listnames)
names(spglist_body) <- listnames
names(spglist_head) <- listnames
 
write(spglist_body$Royal.Female_specific_genes, file="../SSGs/RoyalFemale_Specific_Genes_ThAb.txt", sep="\n")
write(spglist_body$Royal.Male_specific_genes, file="../SSGs/RoyalMale_Specific_Genes_ThAb.txt", sep="\n")
write(spglist_body$Soldier.Female_specific_genes, file="../SSGs/SoldierFemale_Specific_Genes_ThAb.txt", sep="\n")
write(spglist_body$Soldier.Male_specific_genes, file="../SSGs/SoldierMale_Specific_Genes_ThAb.txt", sep="\n")
write(spglist_body$Worker.Female_specific_genes, file="../SSGs/WorkerFemale_Specific_Genes_ThAb.txt", sep="\n")
write(spglist_body$Worker.Male_specific_genes, file="../SSGs/WorkerMale_Specific_Genes_ThAb.txt", sep="\n")
 
write(spglist_head$Royal.Female_specific_genes, file="../SSGs/RoyalFemale_Specific_Genes_Head.txt", sep="\n")
write(spglist_head$Royal.Male_specific_genes, file="../SSGs/RoyalMale_Specific_Genes_Head.txt", sep="\n")
write(spglist_head$Soldier.Female_specific_genes, file="../SSGs/SoldierFemale_Specific_Genes_Head.txt", sep="\n")
write(spglist_head$Soldier.Male_specific_genes, file="../SSGs/SoldierMale_Specific_Genes_Head.txt", sep="\n")
write(spglist_head$Worker.Female_specific_genes, file="../SSGs/WorkerFemale_Specific_Genes_Head.txt", sep="\n")
write(spglist_head$Worker.Male_specific_genes, file="../SSGs/WorkerMale_Specific_Genes_Head.txt", sep="\n")


