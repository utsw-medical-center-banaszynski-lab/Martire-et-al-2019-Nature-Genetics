###MA plot from cuffdiff result####

module load R
#installing packages
install.packages("ggpubr")
library(ggpubr)

#import cuffdiff result
data1 <- read.csv("H33.csv")

#View(data1)
head(data1)

#change sample name
names(data1)[8]<-paste("WT_D0")
names(data1)[9]<-paste("KO_D0")

#extract the columns gene,FPKM,FC,and FDR
data_1 <- data1[,c(3,8,9,10,13)]

#remove FPKM values less than 1 from both samples (optional)
#new_data<-data_1[(data_1$WT_D0>1),]
#new_data2<-new_data[(new_data$KO_D0>1),]

#Remove duplicate genes
no_dupl <- data_1[!duplicated(data_1$gene),]

#Setting up genes column as row names
rownames(no_dupl) <- no_dupl$gene

#select column of interest
no_dupl <- no_dupl[2:5]
#View(no_dupl) #view data

no_dupl_df<- no_dupl[1:2] #selecting sample FPKM for mean
#View(no_dupl_df) #View data

##Calculating row means
mean<-rowMeans(no_dupl_df) #calculating the mean values
#View(mean) #View data
test2<- no_dupl[3:4]
#View(test2)
my_data <- cbind(mean,test2)

#Combining all columns to one using column bind
#my_data<- data.frame(cbind(mean,no_dupl$log2.fold_change.,no_dupl$q_value)) ##doesn't work
#View(my_data) #View data
names<-row.names(my_data) #set rownames
names(my_data) #view row names
#change columnnames
names(my_data)[1]<-paste("baseMean")
names(my_data)[2]<-paste("log2FoldChange")
names(my_data)[3]<-paste("padj")
my_data$log2FoldChange<-as.numeric(as.character(my_data$log2FoldChange)) ##to deal with as.factors error!
#View(my_data)

#plot MA
#Below plot can do magic! Try different types of points, colours, displaying top up & down genes, etc.

ggmaplot(my_data, main = "MA plot", xlab = "Log2 mean of normalized counts", ylab = "Log2FC (KO/WT)", ylim = c(-6, 6),
         fdr = 0.05, fc = 2, size = 1,legend = "top",
         genenames = as.vector(row.names(my_data)),
         top=0,select.top.method = c( "fc"),
         font.label = c("bold", 11),label.rectangle = FALSE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = theme_classic())

dev.off()
