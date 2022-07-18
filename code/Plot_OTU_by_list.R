#Plot OTU by list
library(reshape2)
library(ggplot2)

salivayoung=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Compare_Different_method/Summary/youngsaliva.csv", header = T)
salivaold=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Compare_Different_method/Summary/oldsaliva.csv", header = T)
stoolyoung=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Compare_Different_method/Summary/youngstool.csv", header = T)
stoolold=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Compare_Different_method/Summary/oldstool.csv", header = T)
OTUTABLE1=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/otu_table.csv", header=T, sep=",", row.names=1)
TAXTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/tax.csv", header=T, sep=",", row.names=1)

data=stoolold
path=file.path("~/Desktop/stoolold/")

OTUTABLE=lapply(OTUTABLE1, function(x) x / sum(x))
OTUTABLE=as.data.frame(OTUTABLE)
row.names(OTUTABLE)=row.names(OTUTABLE1)

ti=paste("Sa_Y", c(1:10), sep="_")
t1=paste("Sa_O", c(11:23), sep="_")
t2=paste("Sto_Y", c(1:10), sep="_")
t3=paste("Sto_O", c(11:23), sep="_")
t4=append(ti,t1)
t5=c(t2,t3)
col1=c(t4,t5)
colnames(OTUTABLE) <- col1
Age_Group=append(rep("young",10), rep("old",13))
data <- data[,colSums(is.na(data))<nrow(data)]
group=as.factor(colnames(data))
OTUCount=rowSums(OTUTABLE1)
OTUCount=as.data.frame(OTUCount)

for ( i in group) {
  OTU=data[,grep(i, colnames(data))]
  OTU=OTU[OTU!=""]
  for (j in OTU){
   Target=OTUTABLE[j,24:46]
   tTa=t(Target)
   tTa=melt(tTa)
   tTa=cbind(tTa,Age_Group)
   tax=TAXTABLE[j,]
   otucount=paste("Total_Reads_Counts", OTUCount[j,], sep=":")
   taxnames1=paste(tax$phylum, tax$class,tax$order," \n ", sep="|")
   taxnames2=paste(tax$family,tax$genus,otucount, sep="|")
   taxnames=paste(taxnames1,taxnames2, sep="")
   Bnames=paste(i,j,sep="_")
   Anames=paste(Bnames, taxnames,sep=" \n ")
   mypath=file.path(paste(path,Bnames,".jpeg", sep=""))
   jpeg(filename = mypath, width = 10, height = 8, units = "in", res=800)
   z=ggplot(tTa, aes(x=Var1, y=value, fill=Age_Group)) + geom_bar(stat="identity")+geom_vline(xintercept =10.5, colour="black", linetype="solid") + geom_text(aes(label=sprintf("%.02f %% ", tTa$value*100)), vjust=-0.5, hjust=-0.1, size=2.3, colour = "black", fontface="bold", angle=30) +labs(title=Anames, x="Subject", y = "Reads Count") + theme(text = element_text(size=14),axis.text.x = element_text(angle=90, vjust=1,face = "bold"))
   plot(z)
   dev.off()
   print(Anames)
   }
}




