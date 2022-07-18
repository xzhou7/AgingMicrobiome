#compute different subset in HAM
library("phyloseq")
library("ggplot2")
library("scales")
library("grid")
library("ape")
library("plyr")
OTUTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/otu_table.csv", header=T, sep=",", row.names=1)
TAXTABLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/tax.csv", header=T, sep=",", row.names=1)
SAMPLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/SampleID.csv", header=T, sep=",", row.names=1)
TREE=read_tree("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Sequence Result/OTU95/Physeq_Tree")
#Add taxa data with OTU name
taxmat=as.data.frame(TAXTABLE)
taxmat$species=paste(taxmat$genus, row.names(taxmat), sep="|")
taxmat=as.matrix(taxmat)
#merge data into physeq
physeq=phyloseq((otu_table(as.matrix(OTUTABLE),taxa_are_rows=TRUE)),tax_table(taxmat),sample_data(SAMPLE),TREE)

physeqSaYcore=prune_taxa(SaY, physeq)
physeqSaOcore=prune_taxa(SaO, physeq)
physeqStYcore=prune_taxa(StY, physeq)
physeqStOcore=prune_taxa(StO, physeq)

plot_bar(physeqStOcore, "Xpatientid", "Abundance", "genus", title="Old Stool Core",facet_grid = "phylum~.")



physeqyoungstool=subset_samples(physeq, sampletype=="youngstool")
physeqyoungstool=prune_taxa(taxa_sums(physeqyoungstool)>0, physeqyoungstool)
physeqoldstool=subset_samples(physeq, sampletype=="oldstool")
physeqoldstool=prune_taxa(taxa_sums(physeqoldstool)>0, physeqoldstool)

physeqyoungstoolClostridia=subset_taxa(physeqyoungstool, class=="Clostridia")
physeqoldstoolClostridia=subset_taxa(physeqoldstool, class=="Clostridia")

physeqstoolcore=prune_taxa(StoolOTUlist, physeq)
physeqstoolcoreClostridia=subset_taxa(physeqstoolcore,class=="Clostridia")
physeqstoolcoreClostridia=subset_samples(physeqstoolcoreClostridia, site=="stool")

