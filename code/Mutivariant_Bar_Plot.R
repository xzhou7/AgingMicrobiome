library(RColorBrewer)
sequential <- brewer.pal(12,"Paired")

OTUTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/otu_table.csv", header=T, sep=",", row.names = 1)
otutable=as.matrix(OTUTABLE)
totu=t(OTUTABLE)
totu=as.matrix(totu)

barplot(otutable,names.arg = colnames(otutable), cex.names = 0.7,col = sequential, las=2)
legend("topright", legend = row.names(otutable), fill = sequential, cex=0.5)

salivaOTU=OTUTABLE[SalivaOTUlist,]
dim(salivaOTU)

stoolOTU=OTUTABLE[StoolOTUlist,]
dim(stoolOTU)

SAOTUTABLE=OTUTABLE[,1:23]
SAOTUTABLE=SAOTUTABLE[rowSums(SAOTUTABLE)!=0,]

OLD=intersect(binaryoldsaliva$V1,binaryoldstool$V1)
YOUNG=intersect(binaryyoungsaliva$V1,binaryyoungstool$V1)

STOTUTABLE=OTUTABLE[,24:46]
STOTUTABLE=STOTUTABLE[rowSums(STOTUTABLE)!=0,]
