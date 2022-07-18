OTUTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/otu_table.csv", header=T, sep=",", row.names=1)

OTU=lapply(OTUTABLE, function(x) x/sum(x))
OTU=as.data.frame(OTU)
row.names(OTU)=row.names(OTUTABLE)

#subsample by site
saliva=OTU[,1:23]
stool=OTU[,24:46]


#remove otus not exsited in this site
saliva=saliva[which(!apply(saliva, 1 ,FUN = function(x){all(x == 0)})),]
stool=stool[which(!apply(stool, 1, function(x){all(x==0)})),]

#sum up OTU fraction  
salivasum=cbind(rowSums(saliva[,c(1:10)]),rowSums(saliva[,c(11:23)]))
salivasum=as.data.frame(salivasum)
colnames(salivasum)=c("young", "old")

stoolsum=cbind(rowSums(stool[,c(1:10)]), rowSums(stool[,c(11:23)]))
stoolsum=as.data.frame(stoolsum)
colnames(stoolsum)=c("young", "old")

write.csv(salivasum, "~/Desktop/HAM_salivasum.csv")
write.csv(stoolsum, "~/Desktop/HAM_stoolsum.csv")

stoolsum$combine=stoolsum$young*stoolsum$old
salivasum$combine=salivasum$young*salivasum$old

stoollist=row.names(stoolsum[which(stoolsum$combine==0),])
salivalist=row.names(salivasum[which(salivasum$combine==0),])


library("phyloseq")
library("ggplot2")
library("scales")
library("grid")
library("ape")
library("plyr")

TAXTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/tax.csv", header=T, sep=",", row.names=1)
SAMPLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/SampleID.csv", header=T, sep=",", row.names=1)
TREE=read_tree("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/Physeq_Tree")

#
list=as.character(rfsaliva$V1)
OTUheat=OTUTABLE[list,]
TAXheat=TAXTABLE[list,]

otumat=as.matrix(OTUheat)
taxmat=as.data.frame(TAXheat)
taxmat$species=paste(taxmat$genus, row.names(taxmat), sep="|")
taxmat$stain=row.names(taxmat)
taxmat=as.matrix(taxmat)

#merge data into physeq
OTU=otu_table(otumat,taxa_are_rows=TRUE)
TAX=tax_table(taxmat)
map1=sample_data(SAMPLE)

#random_tree=phy_tree(TREE)
physeq=phyloseq(OTU,TAX,map1)


physeqstool=subset_samples(physeq, site=="stool")
physeqstoolyoung=subset_samples(physeqstool, group== "HY")
physeqrum=subset_taxa(physeqstool, family=="Ruminococcaceae")


p2=plot_heatmap(physeqstool, method="NMDS", sample.label = "sampletype",sample.order="sampletype", taxa.label = "species",
                 low="#FFFFCC", high="#000033", na.value="white",title = "Human Aging Microbiome Stool")
p2 = p2 + theme(axis.text.x = element_text(size = 10, angle = -90, hjust = 0), 
                axis.text.y = element_text(size = 8, angle=0, hjust = 0))
p2

#Network Plot
physeqnet = filter_taxa(physeqstool, function(x) mean(x) > 1e-4, TRUE)
netphyseq=make_network(physeqnet,type="taxa",max.dist = 0.5)
plot_network(netphyseq, physeqnet)
plot_net(physeqrum,type = "taxa", maxdist = 0.3, color="family", point_size=5)


physeqf=transform_sample_counts(physeq, function(x) x/sum(x))
ordu = ordinate(physeq, "PCoA", "bray")
p=plot_ordination(physeqf, ordu, color = "sampletype", shape = "site")+geom_text(mapping = aes(label = patientid), size =8, vjust = 1.5,  face="bold")
p=p+ geom_point(size = 5, alpha = 0.75)
p=p+ scale_colour_brewer(type = "qual", palette = "Set1")
p=p+ ggtitle("PCoA Stool Plot on Bray Curits Distance, HAM")+theme(axis.text=element_text(size=18, colour="black"), axis.title=element_text(size=18, face="bold"), 
axis.title.x=element_text(vjust=-0.4, colour="black"), plot.title=element_text(size=24, face="bold"), 
legend.text=element_text(size=20), legend.title=element_text(size=20, face="bold"))
p

physeqsingleplot <- transform_sample_counts(physeq, function(x) x/sum(x))

for (i in list){
gpsfb = subset_taxa(physeqsingleplot, stain==i)
g=plot_bar(gpsfb,"Xpatientid", title=i)
g
}

for (i in 1:30){
f=as.matrix(OTUheat[i,])
barplot(f, main=row.names(f))
}
