#load required packages
require("groHMM")
library("IRanges")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("GenomicFeatures")
options(mc.cores=getCores(4))
#install.packages("vioplot")
library(vioplot)

######H3K27ac-SI-rep3
WT_rep1 <- readGAlignments("/path-to-nodup-bam/filename1_nodup.bam")
KO_rep1 <- readGAlignments("/path-to-nodup-bam/filename2_nodup.bam")

WT_rep1.gr <- granges(WT_rep1)
KO_rep1.gr <- granges(KO_rep1)

##Assign to variable
WT <- WT_rep1.gr
KO <- KO_rep1.gr

##Combine replicates
##If there are replicates, only then.
#WT <- c(WT_rep1.gr,WT_rep2.gr)
#KO <- c(KO_rep1.gr,KO_rep2.gr)

#take either enhancer or refseq.bed
enhancers <- import("ESC_enhancers.bed", format = "BED")

## Find average
library_WT_129 <- NROW(WT)
library_KO_282 <- NROW(KO)

# Calculate RPKM
#Under enhancers
## Import files RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
rpkm_WT_129 <- countOverlaps(enhancers, WT) / (width(enhancers)/1000 * library_WT_129/1000000)
rpkm_KO_282 <- countOverlaps(enhancers, KO) / (width(enhancers)/1000 * library_KO_282/1000000)
rpkm_enhancers_p300 <- data.frame(rpkm_WT_129,rpkm_KO_282)
head(rpkm_enhancers_p300)
write.table(rpkm_enhancers_p300, file="rpkm_H3K27ac_under_ALLenhancers.xls", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

pdf('boxplot_p300_rep1_rep2_merged_distalDOWN.pdf')
boxplot(rpkm_enhancers_p300, col=(c("blue","red3")), main="p300 merged under distal (DOWN)", ylab="RPKM",outline=FALSE, notch=TRUE, names=c("WT", "H3.3 KO"),yaxt="n", cex.axis=1,las=2,lwd=4,lty=1)
axis(2, at=seq(0,10,.5))

dev.off()

wilcox.test(rpkm_enhancers_p300$rpkm_WT_129, rpkm_enhancers_p300$rpkm_KO_282,conf.int=TRUE)

##End of Script##
