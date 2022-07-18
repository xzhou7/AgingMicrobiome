
##--------------------------------------------------------------------------------------------------------------------##

## load tax file ##

##--------------------------------------------------------------------------------------------------------------------##


taxa   = read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE1/Sequence Results/tax.csv", header=T, sep=",")
otu    = read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE1/Sequence Results/otu_table.csv", header=T, sep=",", row.names=1)
sample = read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE1/Sequence Results/SAMPLE.csv", header=T, sep=",", row.names=1)
  
  
  
osp=lapply(otu, function(x) x/sum(x))
osp=as.data.frame(osp)
osp=cbind(row.names(otu), osp)

taxa$combine = paste(taxa$domain, taxa$phylum, taxa$class, taxa$order, taxa$family, taxa$genus, sep="|")
taxa$combine2=paste(taxa$combine, taxa$X, sep=":")
taxacombine=taxa[,9]
taxacombine=as.matrix(taxacombine)
rownames(taxacombine)=taxa$X
lefsesaliva=merge(osp,taxacombine,by.x="row.names(otu)", by.y="row.names")
row.names(lefsesaliva)=lefsesaliva$V1
data <- lefsesaliva [, -which(colnames(lefsesaliva) %in% c("row.names(otu)" ,"V1"))]
tdata <- t(data)
tdata <- as.data.frame(tdata)

#-----------------------------------------------------------------------------------#

# Subset data for Lefse #

#-----------------------------------------------------------------------------------#

test <- merge(sample, tdata, by="row.names")
ttest <- as.data.frame(t(test))

write.csv(ttest, file="~/Desktop/lefseCRMM.csv")

#-----------------------------------------------------------------------------------#

# Installing a list of Package #

#-----------------------------------------------------------------------------------#
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("agricolae", "gamlss", "glmnet", "inlinedocs", "logging", "optparse", "outliers", "penalized", "psc")



totu=t(lefsesaliva)
totu=as.data.frame(totu)
SAMPLE=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/SampleID.csv", header=T, sep=",", row.names=1)
t_MaAsLin_HAM=merge(SAMPLE, totu, by="row.names")
MaAsLin_HAM=t(t_MaAsLin_HAM)
MaAsLin_HAM=as.data.frame(MaAsLin_HAM)
write.csv(MaAsLin_HAM, "~/Desktop/MaAsLin_HAM.csv")


#lefse for picrust
crmm_module=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE1/Analysis Results/CRMM Picrust/all sample file to work with/Data_out_module_all_more_info_Samples.csv", header=T, sep=",")
crmm_module$combine=paste(crmm_module$function_metabolism, crmm_module$metabolism, crmm_module$HumanN_abundance.Function, sep="|")
row.names(crmm_module)=crmm_module$combine
write.csv(crmm_module, file="~/Desktop/lefseCRMMpicrustallmodule.csv")



crmm_pathway=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE1/Analysis Results/CRMM Picrust/all sample file to work with/Data_out_pathway_all_more_info_Samples.csv", header=T, sep=",")
crmm_pathway$combine=paste(crmm_pathway$function_metabolism, crmm_pathway$metabolism, crmm_pathway$as.character.HumanN_abundance...2.., sep="|")
row.names(crmm_pathway)=crmm_pathway$combine
write.csv(crmm_pathway, file="~/Desktop/lefseCRMMpicrustallpathway.csv")


