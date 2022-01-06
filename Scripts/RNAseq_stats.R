library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork)
library(dplyr)
library(stringr)
library(data.table)
library(tibble)
library(viridisLite)
library(ggproto)
library(Cairo)
library(datasets)
library(PMCMRplus)

ListForStats <- ListByCluster[ListByCluster$features.plot %in% "rna_Grin1",]
ggplot(ListForStats, aes(x=id, y=avg.exp)) + 
  geom_bar(stat = "identity")
ks.test(ListForStats$avg.exp,ListForStats$id,alternative = c("two.sided", "less", "greater"))
#ggplot(ListForStats, aes(x = features.plot)) + geom_histogram(aes(y = avg.exp.scaled)) #+
  #(fun = dnorm, colour = "red",
  #arg = list(mean = mean(airquality$Ozone, na.rm = TRUE),
  #sd = sd(airquality$Ozone, na.rm = TRUE)))
write.csv(sf,"ClusterPoolResults.csv", row.names = TRUE)

#T-test
ListbyClusterAll <- b$data[b$data$id %in% ClusterPoolAll,]
GeneStatResults <- data.frame(t=numeric(), df=numeric(), p.value=numeric(), name=character(), method=character())
GeneStatResults$name <- as.character(GeneStatResults$name)
GeneStatResults$method <- as.character(GeneStatResults$method)
method <- "pct"
for (i in features){
  ListbyGene <- ListbyClusterAll[str_detect(row.names(ListbyClusterAll), i), ]
  ListByCluster1 <- ListbyGene[ListbyGene$id %in% ClusterPool1,]
  ListByCluster2 <- ListbyGene[ListbyGene$id %in% ClusterPool2,]
  if (method == "avg"){
    ClusterStatResults <- t.test(ListByCluster1$avg.exp, ListByCluster2$avg.exp)
  } else if (method == "pct"){
    ClusterStatResults <- t.test(ListByCluster1$pct.exp, ListByCluster2$pct.exp)
  }
  ClusterStatResults <- append(ClusterStatResults, i, after = length(ClusterStatResults))
  x <- ClusterStatResults[[11]]
  NewStatResult <- data.frame(t=ClusterStatResults$statistic, df=ClusterStatResults$parameter, p.value=ClusterStatResults$p.value, name=x, method=method)
  NewStatResult$name <- as.character(NewStatResult$name)
  NewStatResult$method <- as.character(NewStatResult$method)
  GeneStatResults[nrow(GeneStatResults) + 1, ] <- NewStatResult
  #str(ClusterStatResults)
}
GeneStatResults[ , "p.value.adj"] <- as.numeric()
for(i in 1:length(GeneStatResults[,ncol(GeneStatResults)])) {
  GeneStatResults[i ,ncol(GeneStatResults)] <- (GeneStatResults[i,]$p.value)*nrow(GeneStatResults)
}


#Kruskal Walis - The method not chosen
ListByCluster1 <- b$data[b$data$id %in% ClusterPool1,]
ListByCluster2 <- b$data[b$data$id %in% ClusterPool2,]
print(levels(ListByCluster1$features.plot))
WithinStatResults1 <- kruskal.test(avg.exp ~ features.plot, data = ListByCluster1)
WithinStatResults2 <- kruskal.test(avg.exp ~ features.plot, data = ListByCluster1)
print(WithinStatResults1)
print(WithinStatResults2)
PT = pairwise.wilcox.test(ListByCluster1$avg.exp,
                          ListByCluster1$features.plot,
                          p.adjust.method="bonferroni")
print(PT)

#One Way ANOVA - The method chosen and Dunnet T3 Post hoc
ListByCluster1 <- b$data[b$data$id %in% ClusterPool1,]
ListByCluster2 <- b$data[b$data$id %in% ClusterPool2,]
#print(levels(ListByCluster1$features.plot))
WithinStatResults1 <- aov(avg.exp ~ features.plot, data = ListByCluster1)
WithinStatResults2 <- aov(avg.exp ~ features.plot, data = ListByCluster2)
WithinStatResults3 <- aov(pct.exp ~ features.plot, data = ListByCluster1)
WithinStatResults4 <- aov(pct.exp ~ features.plot, data = ListByCluster2)
#print(WithinStatResults1)
#print(WithinStatResults2)
#TukeyHSD(WithinStatResults1)
#TukeyHSD(WithinStatResults2)
PH1 <- dunnettT3Test(WithinStatResults1)$p.value
#g<-t$p.value
PH2 <- dunnettT3Test(WithinStatResults2)$p.value
PH3 <- dunnettT3Test(WithinStatResults3)$p.value
PH4 <- dunnettT3Test(WithinStatResults4)$p.value
write.csv(PH1, file="Excitatory_AVG_PH.csv")
write.csv(PH2, file="Inhibitory_AVG_PH.csv")
write.csv(PH3, file="Excitatory_PCT_PH.csv")
write.csv(PH4, file="Inhibitory_PCT_PH.csv")
summary(WithinStatResults1)
summary(WithinStatResults2)
summary(WithinStatResults3)
summary(WithinStatResults4)
print(t[[1]][["Pr(>F)"]][1])
PH2
PH1 < 0.05

Loop_Cell <- function (x){
  pvalues <- data.frame(p.value=numeric(), x=character(), xend=character())
  pvalues$x <- as.character(pvalues$x)
  pvalues$xend <- as.character(pvalues$xend)
  for(i in 1:ncol(x)){
    for (j in 1:nrow(x)){
      Working_Cell <- x[j,i]
      if (is.na(Working_Cell) == FALSE && Working_Cell <= 0.05){
        wpvalues <- data.frame(p.value=Working_Cell, x=colnames(x)[i], xend=rownames(x)[j])
        wpvalues$x <- as.character(wpvalues$x)
        wpvalues$xend <- as.character(wpvalues$xend)
        pvalues[nrow(pvalues) + 1, ] <- wpvalues
      }
    }
  }
  pvalues <<- pvalues
}
Loop_Cell(g)
pvalues$x[5]
colnames(g)[1]
max(g$rna_Grin2a)
apply(g, 2, max, na.rm=TRUE)

#Function
Run_ANOVA <- function (cluster) {
  Avg_ANOVA <- aov(avg.exp ~ features.plot, data = cluster)
  Avg_Pthoc <- dunnettT3Test(Avg_ANOVA)
  Pct_ANOVA <- aov(pct.exp ~ features.plot, data = cluster)
  Pct_Pthoc <- dunnettT3Test(Pct_ANOVA)
}
Run_ANOVA(ListByCluster1)
Run_ANOVA(ListByCluster2)
#Save the post hoc in a var

#Do the pct.exp post hoc

#Using summary and post hocs (and t-tests above) loop through the data and perform the drawings
#start with the ANOVA, place the pvalue at the side
#Then loop through each not empty cell ignore dashes - and draw based on is it sigificant or no - add funtion for extra stars after
#The line x and x end coordinates will be based on the the row and collumn
#The line y and yend coordinates will be based on a percentage abve the largest of all means. Increase the percentage for each sucessor.

#Multiple Comparisons - Not to be automated
Gene1 <- c("rna_Grin3a")
Gene2 <- c("rna_Grin2b")
#ClusterPoolMC <- ListByCluster2
#MCNumber <- 7
#MCFullResults <- data.frame(t=numeric(), df=numeric(), p.value=numeric(), name=character())
#MCFullResults$name <- as.character(MCFullResults$name)
ListforMC1 <- ClusterPoolMC[ClusterPoolMC$features.plot %in% Gene1,]
ListforMC2 <- ClusterPoolMC[ClusterPoolMC$features.plot %in% Gene2,]
MCResults <- t.test(ListforMC1$pct.exp, ListforMC2$pct.exp)
name=paste(Gene1,Gene2,"pct.exp", "inhibitory", sep = "-")
print(name)
MCFullResults[nrow(MCFullResults) + 1, ] <- data.frame(MCResults$statistic, MCResults$parameter, MCResults$p.value, name)

#Export Table
write.csv(WorkingDataSet, file="MCFullResults.csv")
write.csv(GeneStatResults, file="GeneStatResults l2 corr.csv")

#Adjust P-values
WorkingDataSet <- MCFullResults
WorkingDataSet[ , "p.value.adj"] <- as.numeric()
#print(1:length(WorkingDataSet[,ncol(WorkingDataSet)]))
#(length(WorkingDataSet[["p.value.adj"]]))
#Can also assign p.value.adj to a variable and then use WorkingDataSet["p.value.adj"] to get the length...
for(i in 1:length(WorkingDataSet[,ncol(WorkingDataSet)])) {
  WorkingDataSet[i ,ncol(WorkingDataSet)] <- (WorkingDataSet[i,]$p.value)*MCNumber
}

#Significance Visual
SigVis <- "is.sig"
WorkingDataSet[ , "is.sig"] <- as.character()
SigVisInd <- which(colnames(WorkingDataSet)==SigVis)
#print(length(WorkingDataSet[, SigVisInd]))
for(i in 1:length(WorkingDataSet[, SigVisInd])) {
  if (WorkingDataSet[i,]$p.value.adj < 0.05) {
    WorkingDataSet[i ,SigVisInd] <- "SIG"
  } else {
    WorkingDataSet[i ,SigVisInd] <- "NON"
  }
}

#Find Variance of rna_Grin2a
ListbyGene1 <- ListByCluster1[str_detect(row.names(ListByCluster1), "rna_Grin2a"), ]
ListbyGene2 <- ListByCluster2[str_detect(row.names(ListByCluster2), "rna_Grin2a"), ]
Var1 <- var(ListbyGene1$avg.exp)
Var2 <- var(ListbyGene2$avg.exp)
Var3 <- var(ListbyGene1$pct.exp)
Var4 <- var(ListbyGene2$pct.exp)
Var1
Var2
Var3
Var4
