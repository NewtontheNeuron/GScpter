library(rstudioapi)

#set working directory to the one this file is currently in
setwd(dirname(getActiveDocumentContext()$path))

source("Pre_analysis_functions.R")

# -----------------------------------------------
# Bargraphs for the average expression and percent
# expressed data, and also for both clusterpools

Plot_4_Bar <- function (avg1, pct1, avg2, pct2, c1, c2){
  n_feat <- "rna_Grin2c" #length(features) + 0.1
  
  Bween_pool <- function(method){
    ListbyClusterAll <- b$data[b$data$id %in% ClusterPoolAll,]
    GeneStatResults <- data.frame(t=numeric(), df=numeric(), p.value=numeric(), name=character())
    GeneStatResults$name <- as.character(GeneStatResults$name)
    for (i in features){
      ListbyGene <- ListbyClusterAll[str_detect(row.names(ListbyClusterAll), i), ]
      ListByCluster[[1]] <- ListbyGene[ListbyGene$id %in% ClusterPool1,]
      ListByCluster[[2]] <- ListbyGene[ListbyGene$id %in% ClusterPool2,]
      if (method == "avg"){
        ClusterStatResults <- t.test(ListByCluster[[1]]$avg.exp, ListByCluster[[2]]$avg.exp)
      } else if (method == "pct"){
        ClusterStatResults <- t.test(ListByCluster[[1]]$pct.exp, ListByCluster[[2]]$pct.exp)
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
  avg.lab <- "Average Expression"
  pct.lab <- "%Expressed"
  
  Plot_details <- function (a,b,c,d,e,method){
    ggplot(a, aes(features.plot, d, fill = factor(features.plot))) + 
      geom_col(color="black", show.legend = FALSE) + 
      scale_fill_viridis_d(option = "plasma") + 
      geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width=0.9), width = 0.50) +
      geom_point(data = b, y = c, fill = "white", color="black") +
      scale_y_continuous(expand = c(0,0), limits = c(0, NA)) + 
      labs(x="Gene", y=e) +
      cowplot::theme_cowplot() + 
      scale_x_discrete(labels=clean_label_list) +
      theme(axis.title = element_text(size=20,face="bold")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=20)) +
      theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=20)) +
      annotate("text", label=paste("p < ", Run_ANOVA(b, method, anova.p.val = TRUE)), x = 4.5, y = Position_ANOVA(a, method), size = 8.5, color = "black")
  }
  #Use 1 function to run,correct and draw t tests. Use 5% less of the larget value to plcae the anova pvalue. Export the post hocs as a table.
  Plot_1 <- Plot_details(avg1, c1, exp.def(c1, "avg"), exp.def(avg1, "avg"), avg.lab, "avg") 
  
  Plot_2 <- Plot_details(avg2, c2, exp.def(c2, "avg"), exp.def(avg2, "avg"), avg.lab, "avg")
  Plot_2 <- Plot_Bween(Plot_2, Bween_pool("avg"), "avg", c2)
  
  Plot_3 <- Plot_details(pct1, c1, exp.def(c1, "pct"), exp.def(pct1, "pct"), pct.lab, "pct")
  
  Plot_4 <- Plot_details(pct2, c2, exp.def(c2, "pct"), exp.def(pct2, "pct"), pct.lab, "pct")
  Plot_4 <- Plot_Bween(Plot_4, Bween_pool("pct"), "pct", c2)
  
  Plot <- Plot_1 + Plot_2 + Plot_3 + Plot_4
  return(Plot)
}
Plot <- Plot_4_Bar(AvgExpPar(ListByCluster[[1]]), PctExpPar(ListByCluster[[1]]), AvgExpPar(ListByCluster[[2]]), PctExpPar(ListByCluster[[2]]), ListByCluster[[1]], ListByCluster[[2]])
Plot

save_image('Quad_BarPlot',Plot, width = 5000, height = 4000)
