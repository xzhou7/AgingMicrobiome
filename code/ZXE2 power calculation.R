library(ggplot2)
#library(cowplot)
library(pwr)
library(gtable)
library(grid)
library(reshape2)

cohens_d <- function(x, y) {
  lx <- length(x) - 1
  ly <- length(y) - 1
  md  <- abs(mean(x) - mean(y)) ## mean difference (numerator)
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd) ## common sd computation
  
  cd  <- md/csd ## cohen's d
  print(cd)
  }

cohens_mean <- function(mx, sdx, my, sdy, n1, n2) {
  md <- abs(mx - my)
  n1 <- n1 - 1
  n2 <- n2 - 1
  csd <- (n1 * (sdx^2)) + (n2 * (sdy^2))
  csd <- csd/(n1 + n2)
  csd <- sqrt(csd)
  
  cd <- md/csd
  print(cd)
}

#example of using mean for effect size calculation
cohens_mean(0.011,0.0164,0.0234,0.0164,10,13)

#calculate flu copy limit of detection
x = ((40:55)/10)

#finishied function, x is the limit of detection when sd=0.028 and mean of control group is 5.689
cohen_liofde <- function(x){
  y <- cohens_mean(5.689,0.828,x,0.828,10,10)
  z <- pwr.t.test(d = y, sig.level = 0.05 , power = 0.8, type = c("two.sample"))
  q <- z$n
  q
}

#set x from 4.0 to 5.5
x = ((40:55)/10)
xx = ((400:550)/100)
#plot effect size and flu copy number
plot(xx, cohens_mean(5.689,0.828,xx,0.828,10,13))

for (j in 1:length(x)) 
#create empty datalist to store data from loop
datalist <- data.frame()

#loop x from 4.0 to 5.5, calculate the sample size for each detection
for (i in 1:length(x)) {
  datalist[i,1] <- cohen_lod(x[i]) 
  datalist[i,2] <- x[i]
  datalist[i,3] <- cohens_mean(5.689,0.828,x[i],0.828,10,10)
}

#plot result
p <- ggplot(datalist, aes(x = V1, y = V2)) + geom_point() 
p <- p + geom_smooth(method = "auto", formula = y ~ x,se = TRUE)
p <- p + xlab("sample needed") + ylab("detection of flu viral copy with Power above 80%")
p <- p + geom_hline(yintercept = 5.688856, color = "red")
p1 <- p + geom_hline(yintercept = 4.238104, color = "darkgreen")
p1

p2 <- ggplot(data = datalist, aes(x = V3, y = V2)) + geom_point(color = "red")

ggplot_dual_axis(p1, p2, which.axis = "x")

#example of power calculation from pwr package
pwr.t.test(n=50, sig.level = 0.01 , d = 1, type = c("two.sample"))


plot(pwr.t.test(d = 1.4359903381642514, sig.level = 0.05 , power = 0.8, type = c("two.sample")))

ggplot_dual_axis = function(plot1, plot2, which.axis = "x") {
  
  # Update plot with transparent panel
  plot2 = plot2 + theme(panel.background = element_rect(fill = NA))
  
  grid.newpage()
  
  # Increase right margin if which.axis == "y"
  if(which.axis == "y") plot1 = plot1 + theme(plot.margin = unit(c(0.7, 1.5, 0.4, 0.4), "cm"))
  
  # Extract gtable
  g1 = ggplot_gtable(ggplot_build(plot1))
  
  g2 = ggplot_gtable(ggplot_build(plot2))
  
  # Overlap the panel of the second plot on that of the first
  pp = c(subset(g1$layout, name == "panel", se = t:r))
  
  g = gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)
  
  # Steal axis from second plot and modify
  axis.lab = ifelse(which.axis == "x", "axis-b", "axis-l")
  
  ia = which(g2$layout$name == axis.lab)
  
  ga = g2$grobs[[ia]]
  
  ax = ga$children[[2]]
  
  # Switch position of ticks and labels
  if(which.axis == "x") ax$heights = rev(ax$heights) else ax$widths = rev(ax$widths)
  
  ax$grobs = rev(ax$grobs)
  
  if(which.axis == "x") 
    
    ax$grobs[[2]]$y = ax$grobs[[2]]$y - unit(1, "npc") + unit(0.15, "cm") else
      
      ax$grobs[[1]]$x = ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
  
  # Modify existing row to be tall enough for axis
  if(which.axis == "x") g$heights[[2]] = g$heights[g2$layout[ia,]$t]
  
  # Add new row or column for axis label
  if(which.axis == "x") {
    
    g = gtable_add_grob(g, ax, 2, 4, 2, 4) 
    
    g = gtable_add_rows(g, g2$heights[1], 1)
    
    g = gtable_add_grob(g, g2$grob[[6]], 2, 4, 2, 4)
    
  } else {
    
    g = gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
    
    g = gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b) 
    
    g = gtable_add_grob(g, g2$grob[[7]], pp$t, length(g$widths), pp$b - 1)
    
  }
  
  # Draw it
  grid.draw(g)
  
}


#reproduce figures for HAM 
#set subject number per group
listN <- c(10, 20, 30, 40)
#set significance value
sigN <- c(0.01, 0.001, 0.0001)

datalist <- data.frame()
for (j in 1:length(listN)) {
  print(j)
  for (i in 1:16) {
x <- pwr.t.test(n = (listN[j]), sig.level = 0.05 , d = (3/16*i), type = c("two.sample"))
datalist[i,j] <- x$power
print(i,j)
  }
}
colnames(datalist) <- listN
datamelt <- melt(datalist)
datamelt$delta <- rep(1:16,4) * (3/16)
colnames(datamelt) <- c("Sample.Size.per.group", "Power", "Delta")

p <- ggplot(data = datamelt, mapping = aes(x=Delta, y=Power, color=Sample.Size.per.group)) + geom_point() + geom_line()
p <- p + ggtitle("Two-point Significance=0.05")
p

#Figure 2
#set globle variance
VarN <- c(0.01, 0.0165, 0.05, 0.1, 0.15, 0.2)
#
MeandiffN <- c(0.001,0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.04, 0.06)

x <- cohens_mean(1.02, 0.01, 1.03, 0.01, 10, 13)
datalist3 <- data.frame()

for (j in 1:length(VarN)) {
  print(j)
  for (i in 1:length(MeandiffN)) {
    x <- cohens_mean(1, VarN[j], (1+MeandiffN[i]), VarN[j], 10, 13)
  datalist3[i,j] <-x}
}
colnames(datalist3) <- VarN
datalist3$meandiff <- c(0.001,0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.04, 0.06)
p3list <- paste(1:10, "%", sep="")
datalist3$meandiff<- factor(datalist3$meandiff, levels=c(0.001,0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.04, 0.06))
df3 <- melt(datalist3)
colnames(df3) <- c("Mean.Difference", "Standard.Deviation", "Effect.Size")
p3 <- ggplot(df3, aes(x=Mean.Difference, y=Effect.Size, color=Standard.Deviation)) + geom_point() + geom_line(aes(group=Standard.Deviation))
p3

#figure 1
listnumber0 <- 0:10
listnumber1 <- 12.5 + 10.66667*listnumber0
deltalist <- c(1,2)

datalist2 <- data.frame()
for (i in 1:length(listnumber1)) {
  x <- pwr.t.test(n = (listnumber1[i]/2), sig.level = 0.01 , d = 1, type = c("two.sample"))
  y <- pwr.t.test(n = (listnumber1[i]/3), sig.level = 0.01 , d = 1, type = c("two.sample"))
  datalist2[i,1] <- x$power
  datalist2[i,2] <- y$power
}
colnames(datalist2) <- c("Two.Point", "Three.Point")
data.2.melt <- melt(datalist2)
data.2.melt$N <- rep(round(listnumber1,0), 2)
colnames(data.2.melt) <- c("Design", "Power", "Sample.Size")
p <- ggplot(data = data.2.melt, mapping = aes(x=Sample.Size, y=Power, color=Design)) + xlim(10,125) + ylim(0,1) + geom_point() + geom_line()
p <- p + ggtitle("Delta_1 Significance=0.01")
p
