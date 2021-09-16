library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
levels(clean_neuron_object@active.ident)
features <- c("rna_Cacna1h", "rna_Cacna1g", "rna_Cacna1i")
features <- c("rna_Grin1", "rna_Grin2a", "rna_Grin2b", "rna_Grin2c", "rna_Grin2d", "rna_Grin3a", "rna_Grin3b")
#features <- c("rna_Grin2a", "rna_Grin2b", "rna_Grin2c", "rna_Grin2d")
#features <- c("rna_Grin2b")
a <- DotPlot(clean_neuron_object, features = features) + RotatedAxis()
a
a$data
FeaturePlot(clean_neuron_object, features = features)
cluster.average <- AverageExpression(object=clean_neuron_object, features = features)
head(cluster.average)

DotPlot(object = clean_neuron_object, features = features) + 
  guides(color = guide_colorbar(title = 'Scaled Average Expression')) + 
  theme(axis.text.x = element_text(angle=90))

#Start

b <- DotPlot(clean_neuron_object, features = features)
#b

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

Gene <- b$data$features.plot
Cluster <- b$data$id
AvgExpScaled <- b$data$avg.exp.scaled
markers <- Gene %>% unique()

Plot <- b$data %>% filter(features.plot %in% markers) %>% 
  mutate(`% Expressing` = b$data$pct.exp) %>% 
  ggplot(aes(y=Gene, x = Cluster, color = AvgExpScaled, size = `% Expressing`)) + 
  geom_point() +
  #scale_size(range = c(0, 20)) +
  scale_size_area() +
  scale_color_viridis_c(option = "plasma") + 
  cowplot::theme_cowplot() + 
  #theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1))
Plot

#ClusterPool
ClusterPoolResults <- data.frame(avg.exp=numeric(), pct.exp=numeric(), features.plot=character(), id=character(), avg.exp.scaled=numeric())
ClusterPoolResults$features.plot <- as.character(ClusterPoolResults$features.plot)
ClusterPoolResults$id <- as.character(ClusterPoolResults$id)
# ClusterPoolResults <- ClusterPoolResults %>% add_row(avg.exp=xl[1], pct.exp=xl[2], features.plot=xl[3],id=xl[4], avg.exp.scaled=xl[5]) %>% rownames_to_column(var="23Heaven")
#add_row_now <- function (object, Row_Item){
  #object[nrow(object) + 1, ] = Row_Item
  #object <- object %>% add_row(avg.exp=Row_Item[1], pct.exp=Row_Item[2], features.plot=Row_Item[3],id=Row_Item[4], avg.exp.scaled=Row_Item[5])
  #print(object)
  #row.names(object)[nrow(object)] <- name
#}
#row.names(ClusterPoolResults)[nrow(ClusterPoolResults)] <- "bran"
#ClusterPoolResults <- add_row_now(ClusterPoolResults, xl)
ClusterPool1 <- list("Excit-01", "Excit-02", "Excit-03", "Excit-08", "Excit-09", "Excit-10", "Excit-12", "Excit-14", "Excit-15", "Excit-16", "Excit-18", "Excit-04", "Excit-05", "Excit-13", "Excit-19")
ClusterPool2 <- list("Inhib-01", "Inhib-02", "Inhib-03", "Inhib-04", "Inhib-05", "Inhib-06", "Inhib-07", "Inhib-09", "Inhib-10","Inhib-11","Inhib-12", "Inhib-13")
ClusterPoolAll <- list("Excit-01", "Excit-02", "Excit-03", "Excit-08", "Excit-09", "Excit-10", "Excit-12", "Excit-14", "Excit-15", "Excit-16", "Excit-18", "Inhib-01", "Inhib-02", "Inhib-03", "Inhib-04", "Inhib-05", "Inhib-06", "Inhib-07", "Inhib-09", "Inhib-10","Inhib-11","Inhib-12", "Inhib-13")

#New Cluster Pools
ClusterPool1 <- list("Excit-01", "Excit-02", "Excit-03", "Excit-04", "Excit-05", "Excit-08", "Excit-09", "Excit-10", "Excit-12", "Excit-13", "Excit-14", "Excit-15", "Excit-16", "Excit-18", "Excit-19")
ClusterPool2 <- list("Inhib-01", "Inhib-02", "Inhib-03", "Inhib-04", "Inhib-05", "Inhib-06", "Inhib-07", "Inhib-09", "Inhib-10","Inhib-11","Inhib-12", "Inhib-13")
ClusterPoolAll <- list("Excit-01", "Excit-02", "Excit-03", "Excit-04", "Excit-05", "Excit-08", "Excit-09", "Excit-10", "Excit-12", "Excit-13", "Excit-14", "Excit-15", "Excit-16", "Excit-18", "Excit-19", "Inhib-01", "Inhib-02", "Inhib-03", "Inhib-04", "Inhib-05", "Inhib-06", "Inhib-07", "Inhib-09", "Inhib-10","Inhib-11","Inhib-12", "Inhib-13")

ListByCluster <- b$data[b$data$id %in% ClusterPool1,]
ListByCluster <- b$data[b$data$id %in% ClusterPool2,]
ListbyClusterAll <- b$data[b$data$id %in% ClusterPoolAll,]
#print(ListByCluster[str_detect(row.names(ListByCluster), "rna_Grin1"), ])
for (i in features){
  ListByGene <- ListByCluster[str_detect(row.names(ListByCluster), i), ]
  #print(ListByGene)
  NewListItem <- data.frame(colMeans(ListByGene[1]),colMeans(ListByGene[2]),Col3=i,Col4="Inhibitory",colMeans(ListByGene[5]))
  NewListItem$Col3 <- as.character(NewListItem$Col3)
  NewListItem$Col4 <- as.character(NewListItem$Col4)
  #print(NewListItem)
  #ClusterPoolResults <- ClusterPoolResults %>% add_row_now(ClusterPoolResults, NewListItem)
  ClusterPoolResults[nrow(ClusterPoolResults) + 1, ] <- NewListItem
  #print(ClusterPoolResults)
  row.names(ClusterPoolResults)[nrow(ClusterPoolResults)] <- paste(i, "2")
}
#ListByGene <- ListByCluster[str_detect(row.names(ListByCluster), features[5]), ]
#colMeans(ListByCluster[2])
#Graph
Make_Label <- function(string){
  str_remove(string, "rna_")
}
#Function for pooling the means and re-scaling avg.exp
PoolnShare <- function(cp1, cp2, id1, id2){
  #Defining the new labels
  Make_Label <- function(string){
    str_remove(string, "rna_")
  }
  
  #Defining the full Data Table, it will clear every time
  ClusterPoolResults <- data.frame(avg.exp=numeric(), pct.exp=numeric(), features.plot=character(), id=character(), avg.exp.scaled=numeric(), features.label=character())
  ClusterPoolResults$features.plot <- as.character(ClusterPoolResults$features.plot)
  ClusterPoolResults$id <- as.character(ClusterPoolResults$id)
  ClusterPoolResults$features.label <- as.character(ClusterPoolResults$features.label)
  
  #Separates the data based on the cluster pools
  ListByCluster1 <- b$data[b$data$id %in% cp1,]
  ListByCluster2 <- b$data[b$data$id %in% cp2,]
  
  #For repeating function
  PoolAllRepeat <- 1

  #Average, Scale and save
  PoolAll <- function (cluster, identity){
    for(i in features){
      ListByGene <- cluster[str_detect(row.names(cluster), i), ]
      NewListItem <- data.frame(colMeans(ListByGene[1]), colMeans(ListByGene[2]), Col3=i, Col4=identity, colMeans(ListByGene[5]), Col6=Make_Label(i))
      NewListItem$Col3 <- as.character(NewListItem$Col3)
      NewListItem$Col4 <- as.character(NewListItem$Col4)
      NewListItem$Col6 <- as.character(NewListItem$Col6)
      ClusterPoolResults[nrow(ClusterPoolResults) + 1, ] <<- NewListItem
      if(PoolAllRepeat == 1){
        row.names(ClusterPoolResults)[nrow(ClusterPoolResults)] <<- i
      } else {
        row.names(ClusterPoolResults)[nrow(ClusterPoolResults)] <<- paste(i, PoolAllRepeat)
      }
    }
    PoolAllRepeat <<- PoolAllRepeat + 1
  }
  
  #Running the function PoolAll, (If there were extra pools put them here before the 'Rescale' bellow)
  PoolAll(ListByCluster1, id1)
  PoolAll(ListByCluster2, id2)
  
  #Rescaling average expression based on z scores of the full new dataset
  ClusterPoolResults[, ncol(ClusterPoolResults) + 1] <- 
    data.frame(avg.exp.re.scaled = (ClusterPoolResults[, 1] - colMeans(ClusterPoolResults[1]))/sd(ClusterPoolResults$avg.exp))
  return(ClusterPoolResults)
}
ClusterPoolResults <- PoolnShare(ClusterPool1, ClusterPool2, "Excitatory", "Inhibitory")
#Checking if the clusterpool was generated correctly
#PoolnShare(ClusterPool1, ClusterPool2, "Excitatory", "Inhibitory") == ClusterPoolResults
Gene <- ClusterPoolResults$features.label
Cluster <- ClusterPoolResults$id
AvgExpScaled <- ClusterPoolResults$avg.exp.re.scaled
markers <- Gene %>% unique()

Plot <- ClusterPoolResults %>% filter(features.label %in% markers) %>% 
  mutate(`% Expressing` = ClusterPoolResults$pct.exp) %>% 
  #filter(AvgExpScaled > 1, `% Expressing` > 1) %>% 
  ggplot(aes(y=Gene, x = Cluster, color = AvgExpScaled, size = `% Expressing`)) + #geom_dotplot(dotsize=1.5) + 
  geom_point() + 
  scale_size(range = c(0, 20)) +
  #scale_size_area(max_size = 10) +
  scale_color_viridis_c(option = "plasma") + 
  cowplot::theme_cowplot() + 
  #theme(axis.line  = element_blank()) +
  theme(axis.title = element_text(size=20,face="bold"), legend.key.size=unit(1, "line")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15)) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15)) #+
  #theme(axis.ticks = element_blank()) #+
Plot #+
  #annotate("segment", x = 1.25, xend = 1.75, y = "rna_Grin3b", yend = "rna_Grin3b",
           #colour = "black") +
  #annotate("text", label="*", x = 1.5, y = 7.05, size = 6, color = "black") + 
  #annotate("text", label="+", x = 1.5, y = 6.90, size = 4.4, color = "black")

# Adding Significance Bars
for(i in features) {
  Gene <- i
  return(0)
}

#Saving Immages
ggsave(filename = "Slobo/pooledbargraph.png", plot = Plot,  device = "png", type = "cairo") #, width = 0.1875, height = 27, limitsize = FALSE)
sf <- Plot$data
#scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

  
#Ploting ClusterPoolAll
ListbyClusterAll[, ncol(ListbyClusterAll) + 1] <- data.frame(features.label = Make_Label(ListbyClusterAll$features.plot))
Gene <- ListbyClusterAll$features.label
Cluster <- ListbyClusterAll$id
AvgExpScaled <- ListbyClusterAll$avg.exp.scaled
markers <- Gene %>% unique()

Plot <- ListbyClusterAll %>% filter(features.label %in% markers) %>% 
  mutate(`% Expressing` = ListbyClusterAll$pct.exp) %>% 
  #filter(AvgExpScaled > 1, `% Expressing` > 1) %>% 
  ggplot(aes(y=Gene, x = Cluster, color = AvgExpScaled, size = `% Expressing`)) + #geom_dotplot(dotsize=1.5) + 
  geom_point() + 
  scale_size(range = c(0, 25)) +
  #old raange was to 20
  #scale_size_area(max_size = 10) +
  scale_color_viridis_c(option = "plasma") + 
  cowplot::theme_cowplot() + 
  #theme(text = element_text(family = "serif")) +
  #theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 25), axis.title = element_text(size = 25, face="bold"), legend.key.size = unit(2.1, "line"), legend.text = element_text(size = 17), legend.title = element_text(size = 20)) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 25))
Plot

#>>>>
#Within Cluster Bar graphs
#>>>>
ListByCluster1 <- b$data[b$data$id %in% ClusterPool1,]
ListByCluster2 <- b$data[b$data$id %in% ClusterPool2,]
SE <- function (x) {
  sd(x)/sqrt(length(x))
}
AvgExpPar <- function (x) {
  x %>% group_by(features.plot) %>% summarise(std.err = SE(avg.exp), avg.exp = mean(avg.exp), lower = avg.exp - std.err, upper = avg.exp + std.err)
}
PctExpPar <- function (x) {
  x %>% group_by(features.plot) %>% summarise(std.err = SE(pct.exp), pct.exp = mean(pct.exp), lower = pct.exp - std.err, upper = pct.exp + std.err)
}

Plot_4_Bar <- function (avg1, pct1, avg2, pct2, c1, c2){
  n_feat <- "rna_Grin2c" #length(features) + 0.1
  
  Bween_pool <- function(method){
    ListbyClusterAll <- b$data[b$data$id %in% ClusterPoolAll,]
    GeneStatResults <- data.frame(t=numeric(), df=numeric(), p.value=numeric(), name=character())
    GeneStatResults$name <- as.character(GeneStatResults$name)
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
      NewStatResult <- data.frame(t=ClusterStatResults$statistic, df=ClusterStatResults$parameter, p.value=ClusterStatResults$p.value, name=x)
      NewStatResult$name <- as.character(NewStatResult$name)
      GeneStatResults[nrow(GeneStatResults) + 1, ] <- NewStatResult
    }
    GeneStatResults[ , "p.value.adj"] <- as.numeric()
    for(i in 1:length(GeneStatResults[,ncol(GeneStatResults)])) {
      GeneStatResults[i ,ncol(GeneStatResults)] <- (GeneStatResults[i,]$p.value)*nrow(GeneStatResults)
    }
    return(GeneStatResults)
  }
  #rsource is residing source, i.e. where the significance star resides
  Plot_Bween <- function(wplot, pvalues, method, rsource){
    for(i in 1:nrow(pvalues)){
      Ext_gene <- rsource[rsource$features.plot %in% pvalues$name[i],]
      if(method == "avg"){
        value <- apply(Ext_gene[1], 2, max, na.rm=TRUE) + (apply(Ext_gene[1], 2, max, na.rm=TRUE)*0.10)
      } else if(method == "pct"){
        value <- apply(Ext_gene[2], 2, max, na.rm=TRUE) + (apply(Ext_gene[2], 2, max, na.rm=TRUE)*0.10)
      }
      if(pvalues$p.value[i] < 0.05){
        stars <- "*"
        if(pvalues$p.value[i] < 0.01){
          stars <- "**"
          if(pvalues$p.value[i] < 0.001){
            stars <- "***"
          }
        }
      }
      if(pvalues$p.value.adj[i] < 0.05){
        wplot <- wplot + annotate("text", label = stars, x = pvalues$name[i], y = value, size = 6, colour = "black")
      }
      else {
        #wplot <- wplot + annotate("text", label = "ns", x = pvalues$name[i], y = value, size = 3, colour = "black")
      }
    }
    return(wplot)
  }
  
  Run_ANOVA <- function (cluster, method, anova.p.val=NULL) {
    if(method == "avg"){
      ANOVA <- aov(avg.exp ~ features.plot, data = cluster)
    } else if (method == "pct"){
      ANOVA <- aov(pct.exp ~ features.plot, data = cluster)
    }
    if (isTRUE(anova.p.val)){
      anova <- summary(ANOVA)[[1]][["Pr(>F)"]][1]
      return(signif(anova, digits = 2))
    } else if(is.null(anova.p.val)){
      
    }
    Pthoc <- dunnettT3Test(ANOVA)$p.value
    return(Pthoc)
  }
  
  Loop_Cell <- function (x){
    pvalues <- data.frame(p.value=numeric(), x=character(), xend=character())
    pvalues$x <- as.character(pvalues$x)
    pvalues$xend <- as.character(pvalues$xend)
    for(i in 1:ncol(x)){
      for (j in 1:nrow(x)){
        Working_Cell <- x[j,i]
        if (is.na(Working_Cell) == FALSE && Working_Cell > 0.05){
          wpvalues <- data.frame(p.value=Working_Cell, x=colnames(x)[i], xend=rownames(x)[j])
          wpvalues$x <- as.character(wpvalues$x)
          wpvalues$xend <- as.character(wpvalues$xend)
          pvalues[nrow(pvalues) + 1, ] <- wpvalues
        }
      }
    }
    return(pvalues)
  }
  
  Position_ANOVA <- function(source, method){
    if(method == "avg"){
      pvpos <- apply(source[3], 2, max, na.rm=TRUE) - (apply(source[3], 2, max, na.rm=TRUE)*0.05)
    } else if (method == "pct"){
      increment <- 5
      pvpos <- apply(source[3], 2, max, na.rm=TRUE) - (apply(source[3], 2, max, na.rm=TRUE)*0.05)
    }
    return(pvpos)
  }
  
  highest_mean <- function(source, method, wplot, pvalues){
    if(method == "avg"){
      increment <- 1
      start_point <- apply(source[1], 2, max, na.rm=TRUE) + (apply(source[1], 2, max, na.rm=TRUE)*0.05)
    } else if (method == "pct"){
      increment <- 5
      start_point <- apply(source[2], 2, max, na.rm=TRUE) + (apply(source[2], 2, max, na.rm=TRUE)*0.05)
    }
    for(i in 1:nrow(pvalues)){
      value <- start_point + (i - 1)*increment
      wplot <- wplot + annotate("segment", x = pvalues$x[i], xend = pvalues$xend[i], y = value, yend = value, colour = "black")
    }
    return(wplot)
  }
  
  exp.def <- function (source, method) {
    if(method == "avg"){
      source$avg.exp
    } else if (method == "pct"){
      source$pct.exp
    }
  }
  features_label <- c("Cacna1h", "Cacnag", "Cacnai", "Grin2c", "Grin2d", "Grin3a", "Grin3b")
  avg.lab <- "Average Expression"
  pct.lab <- "%Expressed"
  Plot_details <- function (a,b,c,d,e,method){
    ggplot(a, aes(features.plot, d, fill = factor(features.plot))) + 
      #theme(legend.position = "none") + 
      geom_col(color="black", show.legend = FALSE) + 
      scale_fill_viridis_d(option = "plasma") + 
      #scale_fill_discrete(fill = factor(a$features.plot), guide=FALSE) +
      geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width=0.9), width = 0.50) +
      geom_point(data = b, y = c, fill = "white", color="black") +
      scale_y_continuous(expand = c(0,0), limits = c(0, NA)) + 
      labs(x="Gene", y=e) +
      cowplot::theme_cowplot() + 
      scale_x_discrete(labels=features_label) +
      theme(axis.title = element_text(size=20,face="bold")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=20)) +
      theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=20)) +
      annotate("text", label=paste("p < ", Run_ANOVA(b, method, anova.p.val = TRUE)), x = 4.5, y = Position_ANOVA(a, method), size = 8.5, color = "black")
  }
  #Use 1 function to run,correct and draw t tests. Use 5% less of the larget value to plcae the anova pvalue. Export the post hocs as a table.
  Plot_1 <- Plot_details(avg1, c1, exp.def(c1, "avg"), exp.def(avg1, "avg"), avg.lab, "avg") 
  #Plot_1 <- highest_mean(c1, "avg", Plot_1, Loop_Cell(Run_ANOVA(c1, "avg")))
  Plot_2 <- Plot_details(avg2, c2, exp.def(c2, "avg"), exp.def(avg2, "avg"), avg.lab, "avg") 
  Plot_2 <- Plot_Bween(Plot_2, Bween_pool("avg"), "avg", c2)
  #Plot_2 <- highest_mean(c2, "avg", Plot_2, Loop_Cell(Run_ANOVA(c2, "avg")))
  Plot_3 <- Plot_details(pct1, c1, exp.def(c1, "pct"), exp.def(pct1, "pct"), pct.lab, "pct")
  #Plot_3 <- highest_mean(c1, "pct", Plot_3, Loop_Cell(Run_ANOVA(c1, "pct")))
  Plot_4 <- Plot_details(pct2, c2, exp.def(c2, "pct"), exp.def(pct2, "pct"), pct.lab, "pct")
  Plot_4 <- Plot_Bween(Plot_4, Bween_pool("pct"), "pct", c2)
  #Plot_4 <- highest_mean(c2, "pct", Plot_4, Loop_Cell(Run_ANOVA(c2, "pct")))
  Plot <- Plot_1 + Plot_2 + Plot_3 + Plot_4
  return(Plot)
}
Plot <- Plot_4_Bar(AvgExpPar(ListByCluster1), PctExpPar(ListByCluster1), AvgExpPar(ListByCluster2), PctExpPar(ListByCluster2), ListByCluster1, ListByCluster2)
Plot
t <- AvgExpPar(ListByCluster1)
t
f <- apply(ListByCluster1[2], 2, max, na.rm=TRUE)
f
