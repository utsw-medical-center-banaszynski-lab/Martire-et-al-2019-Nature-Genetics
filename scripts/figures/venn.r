#########Using VennDiagram package in R ##########
#Install the package
install.packages("gdata")
install.packages("VennDiagram")

#Load the package
library(gdata)
require("VennDiagram")
library(gplots)

geneList <- read.xls("GeneList.xlsx", sheet=1, stringsAsFactors=FALSE, header=FALSE)
head(geneList)

# Notice there are empty strings to complete the data frame in column 1 (V1)
tail(geneList)

# To convert this data frame to separate gene lists with the empty strings removed we can use lapply() with our home made  function(x) x[x != ""]
geneLS <- lapply(as.list(geneList), function(x) x[x != ""])

# If this is a bit confusing you can also write a function and then use it in lapply()
removeEMPTYstrings <- function(x) {

 newVectorWOstrings <- x[x != ""]
 return(newVectorWOstrings)

}
geneLS2 <- lapply(as.list(geneList), removeEMPTYstrings)

# You can print the last 6 entries of each vector stored in your list, as follows:
lapply(geneLS, tail)
lapply(geneLS2, tail) # Both methods return the same results

# We can rename our list vectors
names(geneLS) <- c("RNA", "ChIP")

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
#require("VennDiagram")

pdf("venn.pdf")

VENN.LIST <- geneLS
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("RNA", "ChIP"), main="plot title")

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

# To get the list of gene present in each Venn compartment we can use the gplots package
#require("gplots")

a <- venn(VENN.LIST, show.plot=FALSE)

# You can inspect the contents of this object with the str() function
str(a)

# By inspecting the structure of the a object created,
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters <- attr(a,"intersections")

# We can summarize the contents of each venn compartment, as follows:
# in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
lapply(inters, head)
dev.off()
