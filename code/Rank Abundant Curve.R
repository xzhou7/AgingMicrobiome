library(BiodiversityR)
library(vegan)


OTUTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/otu_table.csv", header=T, sep=",", row.names=1)
SAMPLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/SampleID.csv", header=T, sep=",", row.names=1)

saliva=OTUTABLE[,1:23]
stool=OTUTABLE[,24:46]

rsaliva=rrarefy(t(saliva), min(colSums(saliva)))
rowSums(rsaliva)

rstool=rrarefy(t(stool), min(colSums(stool)))
rowSums(rstool)


totuyoungsaliva=rsaliva[1:10,]
totuoldsaliva=rsaliva[11:23,]
totuyoungstool=rstool[1:10,]
totuoldstool=rstool[11:23,]


RankAbun.OSA <- rankabundance(totuoldsaliva)
jpeg("~/Desktop/RAC_Old_Saliva.jpeg",  width = 1080, height = 1080, units = "px", pointsize = 17, quality = 100)
rankabunplot(RankAbun.OSA, scale='proportion', addit=F,xlim=c(0,50), specnames=c(1:10))
dev.off()


RankAbun.OST <- rankabundance(totuoldstool)
jpeg("~/Desktop/RAC_Old_Stool.jpeg",  width = 1080, height = 1080, units = "px", pointsize = 17, quality = 100)
rankabunplot(RankAbun.OST, scale='proportion', addit=F,xlim=c(0,50), specnames=c(1:10))
dev.off()


RankAbun.YSA <- rankabundance(totuyoungsaliva)
jpeg("~/Desktop/RAC_Young_saliva.jpeg",  width = 1080, height = 1080, units = "px", pointsize = 17, quality = 100)
rankabunplot(RankAbun.YSA, scale='proportion', addit=F,xlim=c(0,50), specnames=c(1:10))
dev.off()

RankAbun.YST <- rankabundance(totuyoungstool)
jpeg("~/Desktop/RAC_Young_stool.jpeg",  width = 1080, height = 1080, units = "px", pointsize = 17, quality = 100)
rankabunplot(RankAbun.YST, scale='proportion', addit=F,xlim=c(0,50), specnames=c(1:10))
dev.off()


#B=rankabuncomp(totu, y=SAMPLE, factor='Xpatientid', scale='proportion', legend=F)

for (i in 1:46){
  p=totu[i,]
  RankAbun.1 <- rankabundance(p)
  mypath <- file.path(paste("~/Desktop/",i, sep = ""))
  jpeg(file=mypath)
  rankabunplot(RankAbun.1, scale='proportion', addit=F, specnames=c(1,2,3,4,5),xlim=c(0,+20), main=i)
  dev.off()
}

H=vector("list", length = 46)
for (i in 1:46){
  a<-totu[i,]
  rank<-rankabundance(a)
  b<- row.names(rank)
  H[[i]]<-b
}
H=as.data.frame(H)
SAMPLE$combine=paste(SAMPLE$Xpatientid,SAMPLE$site, sep=".")
colnames(H)=SAMPLE$combine







par(mfcol=c(4,6))


#T test for Renyi Diversity with FDR correction
H=vector("list", length = 8)
qlist=colnames(Renyi.1)
for (i in 1:8) {
  Ttest=t.test(Renyi.1[c(1:10),i],Renyi.1[c(11:23),i])
  Pttest=Ttest$p.value
  H[[i]]=Pttest
}

for (i in 1:8) {
  Ttest=t.test(Renyi.1[c(24:33),i],Renyi.1[c(34:46),i])
  Pttest=Ttest$p.value
  H[[i]]=Pttest
}
H
