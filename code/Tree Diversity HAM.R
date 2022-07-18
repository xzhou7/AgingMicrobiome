#Random Forest Analysis 
#Xin Zhou@Weinstock Lab

#Load Necessary Library for package
library(permute)
library(lattice)
library(vegan)
library(tcltk)
library(BiodiversityR)
library(car)
library(sandwich)
library(RcmdrMisc)
library(splines)
library(Rcmdr)
library(MASS)
library(nlme)
library(mgcv)
library(cluster)
library(RODBC)
library(rpart)
library(effects)
library(mvtnorm)
library(survival)
library(TH.data)
library(multcomp)
library(ellipse)
library(maptree)
library(sp)
library(splancs)
library(spatial)
library(akima)
library(nnet)
library(raster)
library(dismo)
library(rgdal)
library(parallel)
library(gbm)
library(randomForest)
library(foreach)
library(gam)
library(plotmo)
library(plotrix)
library(TeachingDemos)
library(earth)
library(class)
library(mda)
library(kernlab)
library(e1071)
library(sem)
library(rgl)
library(relimp)
library(lmtest)
library(leaps)
library(grid)
library(Formula)
library(ggplot2)
library(Hmisc)
library(colorspace)
library(aplpack)
library(abind)
library(XLConnect)
library(markdown)
library(knitr)



#Rank Abundance: rank the abundance based on t(otumat) 
##RankAbun.1 <- rankabundance(totumat)
##RankAbun.1
##rankabunplot(RankAbun.1, scale="abundance") 
##rankabunplot(RankAbun.1, scale="proportion")

#Input file
OTUTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/otu_table.csv",header=T,sep=",", row.names=1)
SAMPLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/SampleID.csv", header=T, sep=",", row.names=1)



#transpose data to fit package 
tcr=t(OTUTABLE)
tcr=as.matrix(tcr)
#Plot Renyi Diversity
#Renyi.1 = renyiresult(tstoolotu, method="s")
#renyiplot(Renyi.1, evenness=TRUE, col=1)

Renyi.stool = renyicomp(tcr, y=SAMPLE, factor="hypertension", permutations=500, main="HAM Sample Type renyi diversity")


#subset stool sample based on group
tstoolyoung=tstoolotu[1:10,]
tstoolold=tstoolotu[11:23,]

RankAbun.stoolyoung = rankabundance(tstoolyoung)
RankAbun.stoolold = rankabundance(tstoolold)
rankabunplot(RankAbun.stoolold, scale="proportion",labels="old&frail", xlim=c(0,+20),ylim=c(0,+21), col=2, main="Stool Sample Rank Abundance Curve", )
rankabunplot(RankAbun.stoolyoung, scale="proportion", labels="Young", addit=T,col=1)

#subset saliva sample based on group
tsalivayoung=tsalivaotu[1:10,]
tsalivaold=tsalivaotu[11:23,]

RankAbun.salivayoung = rankabundance(tsalivayoung)
RankAbun.salivaold = rankabundance(tsalivaold)
rankabunplot(RankAbun.salivaold, scale="proportion",labels="old&frail", xlim=c(0,+20),ylim=c(0,+13), col=2, main="Saliva Sample Rank Abundance Curve", )
rankabunplot(RankAbun.salivayoung, scale="proportion", labels="Young", addit=T,col=1)
