#Bacterial Evenness Indice
#Xin Zhou 
# Pielou's Evenness Plot

library(vegan)
library(ggplot2)
setwd("~/Dropbox/Weinstock Lab/ZXE/ZXE2/")
getwd()

otu=read.csv("./Sequence Result/OTU95/otu_table.csv", header = T, row.names = 1)
saliva=otu[,1:23]
stool=otu[,24:46]

rsaliva=rrarefy(t(saliva), min(colSums(saliva)))
rstool=rrarefy(t(stool), min(colSums(stool)))


Hstool=diversity(rstool)
Jstool <- Hstool/log(specnumber(rstool))


Hsaliva=diversity(rsaliva)
Jsaliva <- Hsaliva/log(specnumber(rsaliva))

group <- factor(c(rep("young", 10), rep("old", 13)), levels=c("young", "old"))

Evenness=data.frame(group, Jstool, Jsaliva)

write.csv(Evenness, "./Analysis Result/OTU95/RankAbundant/Evenness.csv")

# compare different sample populations across various temperatures
p <- ggplot(EN_saliva, aes(x = group, y = Jsaliva, fill = group)) + geom_boxplot() 
plot (p)

p <- ggplot(EN_stool, aes(x = group, y = Jstool, fill = group)) + geom_boxplot()
plot(p)
p+guides(fill=guide_legend(reverse=TRUE), colour=guide_legend(reverse=TRUE))

plot(Evenness$Jstool, Evenness$Jsaliva, col=group)
