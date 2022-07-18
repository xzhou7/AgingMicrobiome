 
#Load Library
library("phyloseq")
library("ggplot2")
library("scales")
library("grid")
library("ape")
library("plyr")
library("vegan")
#vignette("phyloseq_basics")
#vignette("phyloseq_analysis")

#import file
OTUTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/otu_table.csv", header=T, sep=",", row.names=1)
TAXTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/tax.csv", header=T, sep=",", row.names=1)
SAMPLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/SampleID.csv", header=T, sep=",", row.names=1)
TREE=read_tree("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/Physeq_Tree")
#TREE=read.tree(file = "~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/95otu.tre.txt")
#adjust file based on needs
#OTUTABLE=OTUTABLE[,-grep("rep1", names(OTUTABLE))]
#OTUTABLE=OTUTABLE[,-grep("stool1", names(OTUTABLE))]

otumat=as.matrix(OTUTABLE)
taxmat=as.data.frame(TAXTABLE)
taxmat$species=paste(taxmat$genus, row.names(taxmat), sep="|")
taxmat=as.matrix(taxmat)

#merge data into physeq
OTU=otu_table(otumat,taxa_are_rows=TRUE)
TAX=tax_table(taxmat)
map1=sample_data(SAMPLE)

#random_tree=phy_tree(TREE)
physeq=phyloseq(OTU,TAX,map1,TREE)

#rarefy to even the sequence depth
physeqrare <- rarefy_even_depth(physeq)

#calculate frequency of each taxa
physeqfreq <- transform_sample_counts(physeq,function(x) x/sum(x))

physeqfreqsaliva <- subset_samples(physeqfreq, site=="saliva")
physeqfreqstool <- subset_samples(physeqfreq, site=="stool") 

#Divide by bodysites
physeqraresaliva <- subset_samples(physeqrare, site=="saliva")
physeqrarestool <- subset_samples(physeqrare, site=="stool")

#subset_sample
physeqsaliva=subset_samples(physeq, site=="saliva")
physeqYsaliva=subset_samples(physeqsaliva, sampletype=="youngsaliva")
physeqOsaliva=subset_samples(physeqsaliva, sampletype=="oldsaliva")
physeqstool=subset_samples(physeq,site=="stool")

#Delete Unwanted Taxa
DUT= TAXTABLE["order"]
DUT=subset(DUT, order!="Clostridiales")

DUOTU=rownames(TAXTABLE)
#Setdiff(A, B) means remove all strings in B from A
DUOTU=setdiff(DUOTU, StoolOTUlist)
PhyseqTrimed = prune_taxa(DUOTU, physeq)

listDUT=row.names(DUT)

PhyseqtrimClostridiales = prune_taxa(listDUT, physeq)


#Subset Taxa
physeqOscillibacter=subset_taxa(physeq, genus=="Oscillibacter")
physeqRuminococcus=subset_taxa(physeq, family=="Ruminococcaceae")
physeqAlistipes=subset_taxa(physeq, genus=="Alistipes")
physeqPrevotella=subset_taxa(physeq, genus=="Prevotella")
physeqBacteroides=subset_taxa(physeq, genus=="Bacteroides")
physeqFaecalibacterium=subset_taxa(physeq, genus=="Faecalibacterium")
physeqAlloprevotella=subset_taxa(physeq, genus=="Alloprevotella")
physeqClostridiales=subset_taxa(physeq, order=="Clostridiales")

physeqBacteroidetes=subset_taxa(physeq, phylum=="Bacteroidetes")
physeqBacteroidetesSt=subset_samples(physeqBacteroidetes,site=="stool")

physeqClostridialesST=subset_samples(physeqClostridiales, site=="stool")


physeqClostridialesstool=subset_samples(physeqClostridiales, site=="stool")
PhyseqtrimClostridialesstool=subset_samples(PhyseqtrimClostridiales, site=="stool")

physeqBacteroidetesstool=subset_samples(physeqBacteroidetes, site=="stool")

#Subset_Trimmed_data
physeqTStool=subset_samples(PhyseqTrimed, site=="stool")
# prune OTUs that are not present in at least one sample
physeqTStool = prune_taxa(taxa_sums(physeqTStool) > 0, physeqTStool)

#Heatmap
GP = physeqstool
#Remove OTUs that do not show appear more than 2 times in more than 90% the samples
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 2), A = 0.1 * nsamples(GP))
GP1 = prune_taxa(wh0, GP)
#Remove Taxa that are smaller than 0.05%
GP1 = transform_sample_counts(GP1, function(x) 1e+06 * x/sum(x))
plot_heatmap(GP1,taxa.label="species",sample.order = "sampletype", sample.label = "patientid",distance = "bray")

#Richness
PhyseqRich=prune_taxa(taxa_sums(physeqBacteroidetesSt) > 0, physeqBacteroidetesSt)
plot_richness(stoolcore, x="sampletype", color="age")+geom_boxplot()
estimate_richness(physeqrare, split = TRUE, measures = c("Observed","InvSimpson", "Chao1", "Shannon"))
write.csv(x = estimate_richness(physeq, split = TRUE), file = "~/Desktop/HAM_OTU95_Richness_Plot_rare.csv")

#pick up JSD and Bray as the distance matrix we would like to use
dist_methods <- unlist(distance("list"))
a=dist_methods[c("JSD","vegdist4","betadiver17","dist2","dist3", "vegdist14","vegdist6","vegdist11","betadiver19","betadiver2","betadiver3","betadiver4","betadiver5","betadiver6","betadiver7")]
plist <- vector("list", length(a))
names(plist) = a



for (i in a) {
  # Calculate distance matrix
  iDist <- distance(stoolcore, method = i)
  # Calculate ordination
  iMDS <- ordinate(stoolcore, "PCoA", distance = iDist)
  ## Make plot Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(stoolcore, iMDS, color = "sampletype", shape = "site",label = "patientid")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep = ""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color = sampletype, shape = site, label=patientid))
p = p + geom_point(size = 7)
p = p + geom_text(hjust=0.5, vjust=1.75)
p = p + facet_wrap(~distance, scales = "free")
p = p + ggtitle("Human Aging Microbiome project PCoA Saliva Core Population")
p

#PERMANOVA

#normalize counts by relative abundance
physeq.rel <- transform_sample_counts(physeqfreqstool, function(x) x/sum(x))

df = as(sample_data(physeq.rel), "data.frame")
d = distance(physeq.rel, "bray")
physeq.rel.adonis = adonis(d ~ age, df)
# write.csv(physeq.rel.adonis, "MAM1_PERMANOVA_results_relative abundances.csv")
physeq.rel.adonis
plot(physeq.rel.adonis$aov.tab)

#Ordinations
physeq_ord <- ordinate(physeqfreq, "PCoA", "unifrac", weighted = F)
p4 <- plot_ordination(physeqfreq, physeq_ord, type = "sampletype", color = "sampletype", label="patientid") 
p4 + geom_point(size = 5.5, alpha = 1) + ggtitle("Human Aging Microbiome Unifrac Distance") +  
  theme(axis.text=element_text(size=18, colour="black"), axis.title=element_text(size=18, face="bold"), 
        axis.title.x=element_text(vjust=-0.4, colour="black"), plot.title=element_text(size=24, face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=20, face="bold"))


physeqf=transform_sample_counts(physeq, function(x) x/sum(x))
ordu = ordinate(stoolcore, "PCoA", "bray")
p <- plot_ordination(stoolcore, ordu, color = "sampletype")+scale_colour_manual(values = rhg_cols20) 
p <- p+ scale_shape_manual(values = c(14:19))
#p=p+ geom_text_repel(mapping = aes(label = JCMS.Mouse.ID), size =5,  fontface="bold")
p <- p+ geom_point(size = 7, alpha = 1.0)
#p=p+ scale_colour_brewer(type = "qual", palette = "Set1")
p <- p+ ggtitle("PCoA Plot on Bray_Curtis Distance, MAM")+theme(axis.text=element_text(size=18, colour="black"), axis.title=element_text(size=18, face="bold"), 
                                                                axis.title.x=element_text(vjust=-0.4, colour="black"), plot.title=element_text(size=24, face="bold"), 
                                                                legend.text=element_text(size=20), legend.title=element_text(size=20, face="bold"))
#p <- p + stat_ellipse(mapping=aes(group=Age),type="norm", linetype=2, show.legend=T)+ stat_ellipse(type="t")

p


physeqf=transform_sample_counts(physeqsalivaotuonly, function(x) x/sum(x))
ordu = ordinate(physeqsalivaotuonly, "PCoA", "unifrac", weighted=F)
p=plot_ordination(physeqsalivaotuonly, ordu, color = "sampletype", shape = "site")+geom_text(mapping = aes(label = patientid), size =7, vjust = 1.5,  fontface="bold")
p=p+ geom_point(size = 5, alpha = 0.75)
p=p+ scale_colour_brewer(type = "qual", palette = "Set1")
p=p+ ggtitle("NMDS Stool Plot on Unweighted Unifrac Distance, HAM")+theme(axis.text=element_text(size=18, colour="black"), axis.title=element_text(size=18, face="bold"), 
axis.title.x=element_text(vjust=-0.4, colour="black"), plot.title=element_text(size=24, face="bold"), 
legend.text=element_text(size=20), legend.title=element_text(size=20, face="bold"))
p


physeqnoclostridialesyoungstool=subset_samples(PhyseqtrimClostridiales, sampletype=="youngstool")
physeqnoclostridialesyoungstool=subset_taxa(physeqnoclostridialesyoungstool, order=="Clostridiales")

physeqnoclostridialesoldstool=subset_samples(PhyseqtrimClostridiales,sampletype=="oldstool")
physeqnoclostridialesoldstool=subset_taxa(physeqnoclostridialesoldstool, order=="Clostridiales")

PhyseqTRIMclostridialesoldstool=prune_taxa(binaryoldstool$V1,physeq)

physeqClostridialesyoungstool=subset_samples(physeqClostridiales, sampletype=="youngstool")
physeqClostridialesoldstool=subset_samples(physeqClostridiales, sampletype=="oldstool")

physeqRuminococcusyoungstool=subset_samples(physeqRuminococcus, sampletype=="youngstool")
physeqRuminococcusoldstool=subset_samples(physeqRuminococcus, sampletype=="oldstool")

physeqBsyoungstool=subset_samples(physeqPrevotella, sampletype=="youngsaliva")
physeqBoldstool=subset_samples(physeqPrevotella, sampletype=="oldsaliva")

#Plot Network on taxa
pphyseqncr = filter_taxa(physeqnoclostridialesoldstool, function(x) mean(x) > 1e-6, TRUE)
netphyseq=make_network(pphyseqncr,type="taxa",max.dist = 0.6)
plot_network(netphyseq, pphyseqncr,title = "Stool Old without Binary_Old Clostridiales, Max.Dist=0.6", color="phylum")

plot_net(physeqnoclostridialesoldstool,type = "taxa",maxdist = 0.4, color="phylum")

#Plot network on Sample
ig<-make_network((filter_taxa(physeqrare, function(x) mean(x) > 1e-5, TRUE)), dist.fun="bray", max.dist=0.6)
plot_network(ig, physeq, color="sampletype", label = "patientid", line_weight = 0.5,
             title ="Human Aging Microbiome Network plot/Bray/0.6")



#vennDiagram
library(limma)
physeq100 = subset_samples(physeqrare, patientid=="100")
vennDiagram(object = (vennCounts(otu_table(physeq100))), names = c("saliva", "stool") , 
            lwd = 3, cex = 2, main="Patient 100 VennDiagram")

patientlist=SAMPLE$patientid
patientlist=as.data.frame(patientlist)
patientlist=patientlist[1:23,]


#plot Single Taxa
physeqsingleplot <- transform_sample_counts(physeqsaliva, function(x) x/sum(x))
gpsfb = subset_taxa(physeqsingleplot, phylum=="Firmicutes")
title = "Human Aging Microbiome Project; P/Firmicutes"
plot_bar(gpsfb, "Xpatientid", "Abundance", "genus", title=title,facet_grid = "site~.")


#plot tree

title="Human Aging Microbiome| Stool Clostridia |BootStrap_1000"
plot_tree(physeqstoolClostridia,color="sampletype",size = "abundance",label.tips = "species", text.size = 3,title=title)

#Ordinates On Taxa
GP1= transform_sample_counts(physeqClostridialesST, function(x) 1e+06 * x/sum(x))

GP1=physeqstoolClostridia

phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "family"], sum, na.rm = TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:20]
GP1 = prune_taxa((tax_table(GP1)[, "family"] %in% top5phyla), GP1)

GP.ord <- ordinate(physeq, "NMDS", "bray")
p1 = plot_ordination(GP1, GP.ord, type = "taxa", color = "family")
print(p1)
p1 + facet_wrap(~family)

p3 = plot_ordination(GP1, GP.ord, type = "biplot", color = "class", shape = "sampletype")
# Some stuff to modify the automatic shape scale
#GP1.shape.names = get_taxa_unique(GP1, "phylum")
#GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
#names(GP1.shape) <- GP1.shape.names
#GP1.shape["samples"] <- 16
p3 + facet_wrap(~phylum)


#export distance matrix 
dis=distance(corephyseq, method="bray")
disB=as.matrix(dis)
e <- paste(SAMPLE$Xpatientid, SAMPLE$sampletype, sep=".")
colnames(disB)=rownames(disB)=e
heatmap(log2(disB), Colv = NA, Rowv=NA,symm = T)

disu=UniFrac(corephyseq, weight=T)
disu=as.matrix(disu)
colnames(disu)=rownames(disu)=e
H=log2(disu)
heatmap(log2(disu), Colv=NA, Rowv = NA,symm=T)


disJ=distance(corephyseq, method="jsd")
J=as.matrix(disJ)
colnames(J)=rownames(J)=e
heatmap(log2(J), Colv=NA, Rowv = NA, symm=T)

library(gplots)
LJ=log2(J)
is.na(LJ) <- sapply(LJ, is.infinite)
LJ[is.na(LJ)] <- 0

heatmap.2(LJ,Rowv = FALSE, Colv = FALSE,symkey=T,dendrogram = "none")
write.csv(d, "~/Desktop/HAM_OTU95_Bray.csv")

physeqrare
totu=t(otu_table(physeqrare))

plot_bar(physeqx,"Xpatientid","Abundance", "phylum", title="Saliva Core Microbiome", facet_grid = "site~." )

dis <- disB
df1<-dis[1:10,1:10]
df2<-dis[11:23,11:23]
df3<-dis[1:10, 11:23]
den_d<-density(df1)
den_e<-density(df2)
den_f<-density(df3)

aCol <- rgb(1,0,0,0.2)
bCol <- rgb(0,0,1,0.2)
cCol <- rgb(0,1,0,0.2)

xlim <- range(0.01,den_d$x,den_e$x, den_f$x)
ylim <- range(0,den_d$y, den_e$y,den_f$x)


plot(den_d, xlim = xlim, ylim = ylim, xlab = 'Bray–Curtis dissimilarity Density Plot', main = 'Human Saliva Samples',  panel.first = grid())
polygon(den_d, density = -1, col = aCol)
polygon(den_e, density = -1, col = bCol)
polygon(den_f, density = -1, col = cCol)
legend('topleft',c('Young vs. Young','Aged vs. Aged','Young vs. Aged'), fill = c(aCol, bCol,cCol), bty = 'n',border = 1, cex=1.3)

SAMPLESALIVA <- SAMPLE[1:23,]
attach(SAMPLESALIVA)
dune.dist <- vegdist(t(otu_table(physeqraresaliva)))
dune.ano <- anosim(dune.dist, sampletype, permutations=1000)
summary(dune.ano)



dis <- disB
df1<-dis[24:33,24:33]
df2<-dis[34:46,34:46]
df3<-dis[24:33,34:46]
den_d<-density(df1)
den_e<-density(df2)
den_f<-density(df3)

aCol <- rgb(1,0,0,0.2)
bCol <- rgb(0,0,1,0.2)
cCol <- rgb(0,1,0,0.2)

xlim <- range(0.01,den_d$x,den_e$x, den_f$x)
ylim <- range(0,den_d$y, den_e$y,den_f$x)


plot(den_d, xlim = xlim, ylim = ylim, xlab = 'Bray–Curtis dissimilarity Density Plot', main = 'Human Saliva Samples',  panel.first = grid())
polygon(den_d, density = -1, col = aCol)
polygon(den_e, density = -1, col = bCol)
polygon(den_f, density = -1, col = cCol)
legend('topleft',c('Young vs. Young','Aged vs. Aged','Young vs. Aged'), fill = c(aCol, bCol,cCol), bty = 'n',border = 1, cex=1.3)

SAMPLESTOOL <- SAMPLE[24:46,]
attach(SAMPLESTOOL)
dune.dist <- vegdist(t(otu_table(physeqraresaliva)))
dune.ano <- anosim(dune.dist, sampletype, permutations=1000)
summary(dune.ano)

#Figure S1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#normalize counts by relative abundance
physeq.rel <- transform_sample_counts(physeqfreqstool, function(x) x/sum(x))

df = as(sample_data(physeq.rel), "data.frame")
d = distance(physeq.rel, "bray")
physeq.rel.adonis = adonis(d ~ age, df)
# write.csv(physeq.rel.adonis, "MAM1_PERMANOVA_results_relative abundances.csv")
physeq.rel.adonis
plot(physeq.rel.adonis$aov.tab)

#Ordinations
physeq_ord <- ordinate(physeqfreq, "PCoA", "unifrac", weighted = F)
p4 <- plot_ordination(physeqfreq, physeq_ord, type = "sampletype", color = "sampletype", label="patientid") 
p4 + geom_point(size = 5.5, alpha = 1) + ggtitle("Human Aging Microbiome Unifrac Distance") +  
  theme(axis.text=element_text(size=18, colour="black"), axis.title=element_text(size=18, face="bold"), 
        axis.title.x=element_text(vjust=-0.4, colour="black"), plot.title=element_text(size=24, face="bold"), 
        legend.text=element_text(size=20), legend.title=element_text(size=20, face="bold"))
