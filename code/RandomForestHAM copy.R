#HAM Filter Taxa and Random Forest  
library(randomForest)

OTUTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/otu_table.csv",header=T, sep=",",row.names=1)
t=t(OTUTABLE)

t=as.data.frame(t)
t=lapply(t, function(x) x / sum(x))
t=as.data.frame(t)
row.names(t)=colnames(OTUTABLE)

H=vector("list", length = 548)
for (i in 1:548){
  p=t[,i]
  v=var(p)
  H[[i]]=v
}
H=as.data.frame(H)
H=t(H)
rownames(H)=rownames(OTUTABLE)

Hsub=subset(H, H[,1]>0.005)
Hsub=as.data.frame(Hsub)
Hsub$merge=rownames(Hsub)
OTUTABLE$merge=rownames(OTUTABLE)

Hmerge=merge(Hsub, OTUTABLE, by="merge")
rownames(Hmerge)=Hmerge$merge
Hmerge=Hmerge[,-(1:2)]
r=lapply(Hmerge, function(x) x/sum(x))
r=as.data.frame(r)
row.names(r)=row.names(Hmerge)
r.RF=t(r)
r.RF=as.data.frame(r.RF)
salivarf=r.RF[1:23,]
stoolrf=r.RF[24:46,]

group=c("young","young","young","young","young","young","young","young","young","young","old","old",
        "old","old","old","old","old","old","old","old","old","old","old")
salivarf=cbind(group,salivarf)
stoolrf=cbind(group, stoolrf)

stoolrf=stoolrf[,which(!apply(stoolrf,2,FUN = function(x){all(x == 0)}))]
motustool.rf <- randomForest(group ~ ., data=stoolrf, ntree=1000000, keep.forest=FALSE, importance=TRUE)
varImpPlot(motustool.rf, n.var = 30, main="Stool Sample HAM, ntree=1000000")

salivarf=salivarf[,which(!apply(salivarf,2,FUN = function(x){all(x == 0)}))]
motusaliva.rf <- randomForest(group ~ ., data=salivarf, ntree=1000000, keep.forest=FALSE, importance=TRUE)
varImpPlot(motusaliva.rf, n.var=30, main="Saliva Sample HAM, ntree=1000000")

