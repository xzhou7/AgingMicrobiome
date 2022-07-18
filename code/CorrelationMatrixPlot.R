library(corrplot)

setwd("~/Dropbox/Weinstock Lab/ZXE/ZXE1/")

otu40 <- read.csv("./Analysis Results/Correlationwithcytokines/latefluOTU/top_40_OTU_CRMM.csv", header=T, row.names = 1)
testcytokine <- read.csv("./Lab Data/CRMMImportantSerumdata.csv", header=T, row.names=1)

otu40 <- otu40[!row.names(otu40) %in% "O49",]
cytokine<- testcytokine[!row.names(testcytokine) %in% "O49",]

M <- cor(otu40, cytokine,method = "spearman")


pdf("~/Desktop/Figure4New.pdf", width = 10, height = 12)
corrplot(M, method="circle", tl.pos="lt",tl.col="black",p.mat = res1[[1]],
              sig.level=0.05,insig="blank",tl.cex = 0.8, mar=c(1,1,1,1)) 
dev.off()
?par
#title="Spearman Correlation of most abundant 40 OTUs with Serum Cytokines",

require(FactoMineR) 
# PCA with function PCA

datMy <- cytokine
#read the tab file using the read table function.

pca <- PCA(datMy, scale.unit=TRUE, ncp=5, graph=T)
#scale all the features,  ncp: number of dimensions kept in the results (by default 5)

dimdesc(pca)
#This line of code will sort the variables the most linked to each PC. It is very useful when you have many variables.

cor.mtest <- function(mat1,mat2,conf.level = 0.95){
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n1, n2)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:n1){
    for(j in 1:n2){
      tmp <- cor.test(mat1[,i], mat2[,j], conf.level = conf.level)
      p.mat[i,j]  <- tmp$p.value
      lowCI.mat[i,j]  <- tmp$conf.int[1]
      uppCI.mat[i,j] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(otu40,cytokine,0.95)
res2 <- cor.mtest(testotu,0.99)

corrplot(M, p.mat = res1[[1]], sig.level=0.2)

corrplot.mixed(M)
