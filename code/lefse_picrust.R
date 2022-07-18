

###load tax file
module=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Picrust/all_sample_files_to_work/Data_out_module_all_more_info_Samples.csv", header=T)
pathway=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/Analysis Result/OTU95/Picrust/all_sample_files_to_work/Data_out_pathway_all_more_info_Samples.csv", header=T)
sample=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE2/SampleID.csv", header=T, row.names=1)

module$combine = paste(module$function_metabolism, module$metabolism, module$HumanN_abundance.Function,sep="|")
pathway$combine=paste(pathway$function_metabolism, pathway$metabolism,pathway$as.character.HumanN_abundance...2.., sep="|")

row.names(module)=module$combine
row.names(pathway)=pathway$combine



modulesaliva=module[, 4:26]
modulestool=module[, 27:49]

pathwaysaliva=pathway[,4:26]
pathwaystool=pathway[, 27:49]

e=paste(sample$Xpatientid,sample$sampletype, sep="_")

salivagroup=sample$sampletype[1:23]
stoolgroup=sample$sampletype[24:46]

colnames(modulesaliva)=e[1:23]
colnames(modulestool)=e[24:46]
colnames(pathwaysaliva)=e[1:23]
colnames(pathwaystool)=e[24:46]

modulesaliva=rbind(salivagroup,modulesaliva)
modulestool=rbind(stoolgroup,modulestool)
pathwaystool=rbind(stoolgroup,pathwaystool)
pathwaysaliva=rbind(salivagroup,pathwaysaliva)



write.csv(modulesaliva, file="~/Desktop/HAM_OTU95_modulesaliva.csv")
write.csv(modulestool, file="~/Desktop/HAM_OTU95_modulestool.csv")
write.csv(pathwaysaliva, file="~/Desktop/HAM_OTU95_pathwaysaliva.csv")
write.csv(pathwaystool, file="~/Desktop/HAM_OTU95_pathwaypathwaystool.csv")


taxacombine=taxa[,9]
taxacombine=as.matrix(taxacombine)
rownames(taxacombine)=taxa[,1]
lefsesaliva=merge(osp,taxacombine,by.x="row.names(otu)", by.y="row.names")
lefsesaliva=write.csv(lefsesaliva, file="~/Desktop/lefseCRMMall.csv")

otustool=otu[,71:93]
rownames(otustool)=rownames(otu)

ostp=lapply(otustool, function(x) x/sum(x))
ostp=as.data.frame(ostp)
ostp=cbind(row.names(otustool), ostp)

lefsestool=merge(ostp,taxacombine,by.x="row.names(otustool)", by.y="row.names")
lefsestool=write.csv(lefsestool, file="~/Desktop/lefsestool.csv")

#lefse for picrust
crmm_module=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE1/Analysis Results/CRMM Picrust/all sample file to work with/Data_out_module_all_more_info_Samples.csv", header=T, sep=",")
crmm_module$combine=paste(crmm_module$function_metabolism, crmm_module$metabolism, crmm_module$HumanN_abundance.Function, sep="|")
row.names(crmm_module)=crmm_module$combine
write.csv(crmm_module, file="~/Desktop/lefseCRMMpicrustallmodule.csv")



crmm_pathway=read.csv("~/Dropbox/Weinstock Lab/ZXE/ZXE1/Analysis Results/CRMM Picrust/all sample file to work with/Data_out_pathway_all_more_info_Samples.csv", header=T, sep=",")
crmm_pathway$combine=paste(crmm_pathway$function_metabolism, crmm_pathway$metabolism, crmm_pathway$as.character.HumanN_abundance...2.., sep="|")
row.names(crmm_pathway)=crmm_pathway$combine
write.csv(crmm_pathway, file="~/Desktop/lefseCRMMpicrustallpathway.csv")


