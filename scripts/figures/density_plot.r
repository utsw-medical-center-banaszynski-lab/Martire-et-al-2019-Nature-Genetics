##Density plots in R.

####DensityPlots
##Comparing two columns of data.

##For RPKM FC of WT/KO under Enhancers and Promoters.
Data1<-read.csv("/path-to/yourfolder/WT.csv")
Data2<-read.csv("/path-to/yourfolder/KO.csv")
head(Data1)
head(Data2)
pdf("densityplot_FC_RPKM_at_enh_prom.pdf")
plot(density(Data1$WT_FC),col="darkorchid4", lty=1, lwd=4, main="FC of WT/KO H3K27ac at enh & prom",xlab="RPKM", ylim=c(0,1), xlim=c(-5,5))
lines(density(Data2$KO_FC),col="forestgreen", lty=1, lwd=4)
legend('topright',c('enhancer','promoter'), fill = c("darkorchid4","forestgreen"), bty = 'n', border = NA)
dev.off()
##To test significance & get a statistic: Two-sample Kolmogorov-Smirnov test
ks.test(Data1$WT_FC,Data2$KO_FC)

##For comparing fold-change difference in FPKM between WT (D4/D0) & KO (D4/D0) at NNG.
Data1<-read.csv("/path-to/yourfolder/WT.csv")
Data2<-read.csv("/path-to/yourfolder/KO.csv")
head(Data1)
head(Data2)
pdf("densityplot_FC_FPKM_at_NNG.pdf")
plot(density(Data1$WT_FC),col="darkorchid4", lty=1, lwd=4, main="FC difference in FPKM in WT & KO D4/D0",xlab="RPKM", ylim=c(0,0.70), xlim=c(-10,15))
lines(density(Data2$KO_FC),col="darkorange", lty=1, lwd=4)
legend('topright',c('WT(D4/D0)','KO(D4/D0)'), fill = c("darkorchid4","darkorange"), bty = 'n', border = NA)
dev.off()
##To test significance & get a statistic: Two-sample Kolmogorov-Smirnov test
ks.test(Data1$WT_FC,Data2$KO_FC)