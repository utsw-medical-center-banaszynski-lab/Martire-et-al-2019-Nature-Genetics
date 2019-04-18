# Need these libraries
library(reshape2)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(ggpubr)


tab <- read.table("rpkm.xls")
df <- melt(tab)
# Log scale the RPKM value. + 1 keeps the 0 as 0.
df$log <- log2(df$value + 1)

boxes <- ggplot(df, aes(variable, log,fill=variable)) +
			geom_boxplot(notch =TRUE)+
			theme_classic()+
			scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
			geom_signif(comparisons = list(c("rpkm_KO", "rpkm_WT")),map_signif_level = TRUE,y_position=16)+
			labs(title="WT vs KO",x="Treatment", y = "log2(RPKM)")+
			guides(fill=FALSE)+
			theme(axis.text.x = element_text(angle = 45, hjust = 1))

violins <- ggplot(df, aes(variable, log,fill=variable)) +
			geom_violin(trim=FALSE)+
			geom_boxplot(width=0.1)+
			theme_classic()+
			scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
			geom_signif(comparisons = list(c("rpkm_KO", "rpkm_WT")),map_signif_level = TRUE,y_position=16)+
			labs(title="WT vs KO",x="Treatment", y = "log2(RPKM)")+
			guides(fill=FALSE)+
			theme(axis.text.x = element_text(angle = 45, hjust = 1))

wilcox.test(tab$rpkm_WT, tab$rpkm_KO,conf.int=TRUE)

plot <- plot_grid(boxes, violins,labels=c("A", "B"), ncol = 2,align = "v")
save_plot("Vioplot.pdf", plot, ncol = 2,base_height=6,base_width=3.5)
