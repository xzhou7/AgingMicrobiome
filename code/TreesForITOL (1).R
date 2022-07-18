library(vegan)
library(mgcv)
library(nlme)
library(MASS)
library(labdsv)
library(amap)
library(ctc)

InputLocation="~/Dropbox/Weinstock Lab/ZXE/ZXE8/Sequence Results/"
FileName = "otu_fraction.csv"
ForInput = paste(InputLocation, FileName, sep="")

#make tree
hmpdat=read.csv(file=ForInput, header=TRUE, sep = ",", quote = "\"", dec = ".", stringsAsFactors=FALSE)
row.names(hmpdat) <- hmpdat$genus
hmpdat$X.1 <- hmpdat$X <- hmpdat$genus <- NULL

test=colSums(hmpdat)
hmpdat=hmpdat[,test>0]
hmpdat=prop.table(as.matrix(hmpdat),2)*100
thmpdat=t(hmpdat)
dfhmpdat<-data.frame(thmpdat)
hmpdist<-dsvdis(thmpdat,index="bray/curtis")
hmpclust<-hclust(hmpdist, method="average"  )
write(hc2Newick(hmpclust),file="Genus.class.manhattan.aver.newick")

#make annotation file
hmpdat=read.csv(file=ForInput, header=TRUE, sep = ",", quote = "\"", dec = ".", stringsAsFactors=FALSE)
row.names(hmpdat) <- hmpdat$genus
hmpdat$X.1 <- hmpdat$X <- hmpdat$genus <- NULL

hmpdat=prop.table(as.matrix(hmpdat),2)*100
thmpdat=t(hmpdat)
cs=colSums(thmpdat)
cs.d=data.frame(cs)
cs.d=t(cs.d)
names(cs.d)=names(thmpdat)
total=rbind(thmpdat,cs.d)
n=nrow(total)
o=order(total[n,],decreasing=TRUE)
re.o=total[,o]
re.o=re.o[-n,]
ifelse(ncol(re.o) >= 25, (final=re.o[,c(1:25)]), (final=re.o))
# final=re.o[,c(1:25)]
COLORS=c("#00FFFF","#000000","#0000FF","#FF00FF","#778899","#008000","#800000","#32CD32","#87CEFA","#808000","#800080","#ff0000","#C0C0C0","#008080","#FF6347","#FFFF00","#7CFC00","#000080","#FFD700","#1E90FF","#F778A1","#A0522D","#E0B0FF","#493D26","#FFA62F")
c.final=rbind(COLORS,final)
write.table(c.final,"Genus.dataset.txt",sep="\t",col.names=NA)

#Create color label table
ColorLabel <- data.frame(Var1 = colnames(hmpdat), Label = "label", Color = "black", Notes = "")
ColorLabel$Color <- as.character(ColorLabel$Color)
ColorLabel$Color[grepl("*3..$", ColorLabel$Var1)] <- "#EC2F34"
ColorLabel$Color[grepl("*2..$", ColorLabel$Var1)] <- "#09630B"
ColorLabel$Color[grepl("*1..$", ColorLabel$Var1)] <- "030F89"
write.table(ColorLabel, "ColorLabel.txt", sep = "\t", col.names=NA)

colorlabel=read.csv("sample.csv",header=TRUE,row.names=1,stringsAsFactors=FALSE)
colorlabel1=colorlabel$walnut
colorlabel1<- replace(colorlabel1, which(colorlabel1 == "0walnut"), "#FF7F50")
colorlabel1<- replace(colorlabel1, which(colorlabel1 == "7walnut"), "#7FFF00")

colorlabel2=colorlabel$cancer
colorlabel2<- replace(colorlabel2, which(colorlabel2 == "AOM"), "#EC2F34")
colorlabel2<- replace(colorlabel2, which(colorlabel2 == "NaCl"), "#09630B")

color=cbind(row.names(colorlabel), colorlabel1,colorlabel2)
write.table(color, "ColorLabel.txt", sep = "\t", col.names=NA)

