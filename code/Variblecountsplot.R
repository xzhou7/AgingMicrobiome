

otuSaO=otu_table(physeqSaO_Saliva)
otuSaO=as.data.frame(otuSaO)

MatVar <- function(x, dim = 1) {
  if(dim == 1){
    rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
  } else if (dim == 2) {
    rowSums((t(x) - colMeans(x))^2)/(dim(x)[1] - 1)
  } else stop("Please enter valid dimention")
}
Varible <- MatVar(otuSaO, 1)
Counts <- rowSums(otuSaO!=0)

otuSaOplot=cbind(Varible,Counts)
otuSaOplot <- as.data.frame(otuSaOplot)
library(ggrepel)
p <- ggplot (otuSaOplot, aes(Varible,Counts))
p = p + geom_point()
p = p + geom_text_repel(aes(label=row.names(otuSaOplot),size=2))
p
p<- NULL
