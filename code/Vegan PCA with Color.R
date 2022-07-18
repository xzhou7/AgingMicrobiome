library(vegan)
data(dune)
data(dune.env)
dune.ca <- cca(dune)

colvec <- c("red", "green", "mediumblue","yellow")

ef <- envfit(dune.ca, dune.env, permutations = 999)
plot(dune.ca, display = "sites", type = "p")
with(dune.env, ordiellipse(dune.ca, Management, kind = "se", conf = 0.95))
with(dune.env, ordispider(dune.ca, Management, col = "blue", label= TRUE))
with(dune.env, ordihull(dune.ca, Management, col="blue", lty=2))


dune.otutable=t(OTUTABLE)
dune.otutable=as.data.frame(dune.otutable)
dune.stool=dune.otutable[24:46,]
dune.saliva=dune.otutable[1:23,]
SAMPLE.stool=subset(SAMPLE, SAMPLE$site=="stool")
SAMPLE.saliva=subset(SAMPLE, SAMPLE$site=="saliva")

dune.ca.saliva<- cca(dune.saliva)
ef.otu <- envfit(dune.ca.otu, SAMPLE, permutations = 999)
plot(dune.saliva.rda, display = "sites", type = "p")
with(SAMPLE, ordiellipse(dune.ca.otu, group, kind = "se", conf = 0.95))
with(SAMPLE, ordispider(dune.ca.otu, group, col = "blue", label= TRUE))
with(SAMPLE, ordihull(dune.ca.otu, group, col="blue", lty=2))

dune.saliva.rda=rda(dune.saliva)
ef.saliva.otu <- envfit(dune.saliva.rda, SAMPLE.saliva, permutations = 999)
plot(dune.saliva.rda, display = "sites", type = "p")
with(SAMPLE.saliva, ordiellipse(dune.saliva.rda, group, kind = "se", conf = 0.95))
with(SAMPLE.saliva, ordispider(dune.saliva.rda, group, col = "blue", label= TRUE))
with(SAMPLE.saliva, ordihull(dune.saliva.rda, group, col="red", lty=2))



dune.stool.rda=rda(dune.stool)
plot(dune.stool.rda, display = "sites", type = "p")
with(SAMPLE.stool, ordiellipse(dune.stool.rda, group, kind = "se", conf = 0.95))
with(SAMPLE.stool, ordispider(dune.stool.rda, group, col = colvec, label= TRUE))
with(SAMPLE.stool, ordihull(dune.stool.rda, group, lty=2))



