#Compare Different Method

lefsesaliva=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Compare_Different_method/lefsesaliva.csv",header=F)
lefsestool=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Compare_Different_method/lefsestool.csv",header=F)
rfsaliva=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Compare_Different_method/rfsaliva.csv",header=F)
rfstool=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Compare_Different_method/rfstool.csv",header=F)
binarysaliva=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Compare_Different_method/binarysaliva.csv",header=F)
binarystool=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Compare_Different_method/binarystool.csv",header=F)

lefseyoungsaliva=lefsesaliva[lefsesaliva$V2=="LefseYoungSaliva",]
lefseoldsaliva=lefsesaliva[lefsesaliva$V2=="LefseOldSaliva",]

lefseyoungstool=lefsestool[lefsestool$V2=="LefseYoungStool",]
lefseoldstool=lefsestool[lefsestool$V2=="lefseoldstool",]

rfyoungsaliva=rfsaliva[rfsaliva$V2=="RFYoungSaliva",]
rfoldsaliva=rfsaliva[rfsaliva$V2=="RFOldSaliva",]

rfyoungstool=rfstool[rfstool$V2=="RFYoungStool",]
rfoldstool=rfstool[rfstool$V2=="RFOidStool",]

binaryyoungsaliva=binarysaliva[binarysaliva$V2=="BinarySalivaYoung",]
binaryoldsaliva=binarysaliva[binarysaliva$V2=="BinarySalivaOld",]

binaryyoungstool=binarystool[binarystool$V2=="BinaryStoolYoung",]
binaryoldstool=binarystool[binarystool$V2=="BinaryStoolOld",]

salivaoldlist <- merge((merge(lefseoldsaliva, rfoldsaliva, by="V1", all = T)), binaryoldsaliva, by="V1", all=T)
stoololdlist <- merge((merge(lefseoldstool, rfoldstool, by="V1", all = T)), binaryoldstool, by="V1", all=T)

salivayounglist <- merge((merge(lefseyoungsaliva, rfyoungsaliva, by="V1", all = T)), binaryyoungsaliva, by="V1", all=T)
stoolyounglist <- merge((merge(lefseyoungstool, rfyoungstool, by="V1", all = T)), binaryyoungstool, by="V1", all=T)


TAXTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/tax.csv", header=T, sep=",", row.names=1)
TAXTABLE$V1=row.names(TAXTABLE)
salivayoungtax <- merge (salivayounglist, TAXTABLE, by="V1")
salivaoldtax  <- merge (salivaoldlist, TAXTABLE, by="V1")

stoolyoungtax <- merge(stoolyounglist, TAXTABLE, by="V1")
stoololdtax <- merge(stoololdlist, TAXTABLE, by="V1")



write.csv(salivayoungtax, "~/Desktop/SalivaYounglistcombinedtax.csv")
write.csv(salivaoldtax, "~/Desktop/SalivaOldlistcombinedtax.csv")
write.csv(stoolyoungtax,"~/Desktop/StoolYounglistcombinedtax.csv")
write.csv(stoololdtax,"~/Desktop/StoolOldlistcombinedtax.csv")

SaY=as.character(salivayounglist$V1)
SaO=as.character(salivaoldlist$V1)
SalivaOTUlist=c(SaY,SaO)

StY=as.character(stoolyounglist$V1)
StO=as.character(stoololdlist$V1)
StoolOTUlist=c(StY, StO)

physeqSaY=prune_taxa(SaY, physeq)
physeqSaY_Saliva=subset_samples(physeqSaY, site=="saliva")
sum(rowSums(otu_table(physeqSaY_Saliva)))

physeqSaO=prune_taxa(SaO, physeq)
physeqSaO_Saliva=subset_samples(physeqSaO, site=="saliva")
sum(rowSums(otu_table(physeqSaO_Saliva)))

physeqStY=prune_taxa(StY, physeq)
physeqStY_Stool=subset_samples(physeqStY, site=="stool")
sum(rowSums(otu_table(physeqStY_Stool)))

physeqStO=prune_taxa(StO, physeq)
physeqStO_stool=subset_samples(physeqStO, site=="stool")
sum(rowSums(otu_table(physeqStO_stool)))

physeqBinaryYoungSaliva=prune_taxa(as.character(binaryyoungsaliva$V1), physeq)
physeqBinaryYoungSaliva_Saliva=subset_samples(physeqBinaryYoungSaliva, site=="saliva")
sum(rowSums(otu_table(physeqBinaryYoungSaliva_Saliva)))

physeqBinaryOldSaliva=prune_taxa(as.character(binaryoldsaliva$V1), physeq)
physeqBinaryOldSaliva_Saliva=subset_samples(physeqBinaryOldSaliva, site=="saliva")
sum(rowSums(otu_table(physeqBinaryOldSaliva_Saliva)))


physeqBinaryYoungStool=prune_taxa(as.character(binaryyoungstool$V1), physeq)
physeqBinaryYoungStool_Stool=subset_samples(physeqBinaryYoungStool, site=="stool")
sum(rowSums(otu_table(physeqBinaryYoungStool_Stool)))

physeqBinaryOldStool=prune_taxa(as.character(binaryoldstool$V1), physeq)
physeqBinaryOldStool_Stool=subset_samples(physeqBinaryOldStool, site=="stool")
sum(rowSums(otu_table(physeqBinaryOldStool_Stool)))

physeqrfYoungSaliva=prune_taxa(as.character(rfyoungsaliva$V1), physeq)
physeqrfYoungSaliva_Saliva=subset_samples(physeqrfYoungSaliva, site=="saliva")
sum(rowSums(otu_table(physeqrfYoungSaliva_Saliva)))

physeqrfOldSaliva=prune_taxa(as.character(rfoldsaliva$V1), physeq)
physeqrfOldSaliva_Saliva=subset_samples(physeqrfOldSaliva, site=="saliva")
sum(rowSums(otu_table(physeqrfOldSaliva_Saliva)))


physeqrfYoungStool=prune_taxa(as.character(rfyoungstool$V1), physeq)
physeqrfYoungStool_Stool=subset_samples(physeqrfYoungStool, site=="stool")
sum(rowSums(otu_table(physeqrfYoungStool_Stool)))

physeqrfOldStool=prune_taxa(as.character(rfoldstool$V1), physeq)
physeqrfOldStool_Stool=subset_samples(physeqrfOldStool, site=="stool")
sum(rowSums(otu_table(physeqrfOldStool_Stool)))


physeqlefseYoungSaliva=prune_taxa(as.character(lefseyoungsaliva$V1), physeq)
physeqlefseYoungSaliva_Saliva=subset_samples(physeqlefseYoungSaliva, site=="saliva")
physeqlefseYoungSaliva_Saliva
sum(rowSums(otu_table(physeqlefseYoungSaliva_Saliva)))

physeqlefseOldSaliva=prune_taxa(as.character(lefseoldsaliva$V1), physeq)
physeqlefseOldSaliva_Saliva=subset_samples(physeqlefseOldSaliva, site=="saliva")
physeqlefseOldSaliva_Saliva
sum(rowSums(otu_table(physeqlefseOldSaliva_Saliva)))


physeqlefseYoungStool=prune_taxa(as.character(lefseyoungstool$V1), physeq)
physeqlefseYoungStool_Stool=subset_samples(physeqlefseYoungStool, site=="stool")
physeqlefseYoungStool_Stool
sum(rowSums(otu_table(physeqlefseYoungStool_Stool)))

physeqlefseOldStool=prune_taxa(as.character(lefseoldstool$V1), physeq)
physeqlefseOldStool_Stool=subset_samples(physeqlefseOldStool, site=="stool")
physeqlefseOldStool_Stool
sum(rowSums(otu_table(physeqlefseOldStool_Stool)))


OTUTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/otu_table.csv",header=T, sep=",",row.names=1)
list=as.character(lefseyoungstool$V1)
OTUheat=OTUTABLE[list,]

par(mfrow=c(2, 1))

OTUheat=as.matrix(OTUheat)
a=c(c(100:109), c(200:209), c(300:302))
b=paste("saliva",a,sep="_")
c=paste("stool", a, sep="_")
d=c(b, c)

par(mfrow=c(3, 2))
for (i in 1:6){
  f=as.matrix(OTUheat[i,])
  barplot(f, main=row.names(f))
}


