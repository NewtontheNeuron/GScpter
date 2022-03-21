Bween_pool <- function(method, c1, c2){

    GeneStatResults <- data.frame(t=numeric(), df=numeric(), p.value=numeric(), name=character())
    GeneStatResults$name <- as.character(GeneStatResults$name)

# Start a counter for the total number of comparisons
tcomparisons <- 0
# Function computes a t-test between two cluster pools
# Takes the method = avg or pct, and list by cluster 1 and 2
Bween_pool <- function(method, c1, c2){
    # Use all_cell_roster and cell_roster inherited from pre_analysis_functions
    GeneStatResults <- data.frame(t = numeric(), df = numeric(), p.value = numeric(), features.label = character())
    GeneStatResults$features.label <- as.character(GeneStatResults$features.label)

    for (i in features){ 

      if (method == "avg"){
        ClusterStatResults <- t.test(c1$avg.exp, c2$avg.exp)
      } else if (method == "pct"){
        ClusterStatResults <- t.test(c1$pct.exp, c2$pct.exp)
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

#run anova to get p value 
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

#create individual plot with wanted options.
# Re working plot details to work with CPR_new
Plot_details <- function (pro_data, clusterpool, clusterpool_exp, method_exp, title, method, label, y_lim){
  plot <- ggplot(pro_data, aes(features.label, method_exp, fill = factor(features.label))) + 
    geom_col(color = "black", show.legend = FALSE) + 
    scale_fill_viridis_d(option = "plasma") + 
    #geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width=0.9), width = 0.50) +
    geom_point(data = clusterpool, y = clusterpool_exp, fill = "white", color="black") +
    # Removing the points.
    scale_y_continuous(expand = c(0,0), limits = c(0, y_lim)) + 
    labs(x = "Gene", y = title) +
    cowplot::theme_cowplot() + 
    scale_x_discrete(labels = features_no_key) +
    theme(axis.title = element_text(size = 20, face = "bold")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 20)) +
    theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 20)) +
    annotate("text", label=paste("p < ", Run_ANOVA(clusterpool, method, anova.p.val = TRUE)), 
             x = 4.5, y = Position_ANOVA(pro_data, method), size = 8.5, color = "black")
  
  #add breaks if needeed
  #if (FALSE & method = "avg"){
  #  plot <- plot + 
  #    scale_y_continuous(breaks = (seq(0, 100, by = 20)))
  #}
  
    clean_label_list <- returnCleanLabelList()

    plot <- ggplot(avg, aes(features.plot, method_exp, fill = factor(features.plot))) + 
            geom_col(color="black", show.legend = FALSE) + 
            scale_fill_viridis_d(option = "plasma") + 
            geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width=0.9), width = 0.50) +
            scale_y_continuous(expand = c(0,0), limits = c(0, y_lim)) + 
            labs(x="Gene", y=title) +
            cowplot::theme_cowplot() + 
            scale_x_discrete(labels=clean_label_list) +
              theme(axis.title = element_text(size=20,face="bold")) +
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=20)) +
              theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size=20)) +
            annotate("text", label=paste("p < ", Run_ANOVA(clusterpool, method, anova.p.val = TRUE)), 
              x = 4.5, y = Position_ANOVA(avg, method), size = 8.5, color = "black")
    
    #add a category title for comparisons between the two clusterpools.
    if (method == "avg"){
      plot <- plot + ggtitle(label) +
        theme(plot.title = element_text(hjust = 0.5, size = 34))
    }
    
    return(plot)
    
}

# Bargraphs for the average expression and percent expressed data, and also for both clusterpools
Plot_4_Bar <- function (avg1, pct1, avg2, pct2, c1, c2, labels){

  Plot_1 <- Plot_details(avg1, c1, c1$avg.exp, avg1$avg.exp, "Average Expression", "avg", labels[1], 12) 

  Plot_2 <- Plot_details(avg2, c2, c2$avg.exp, avg2$avg.exp, "Average Expression", "avg", labels[2], 12)
  Plot_2 <- Plot_Bween(Plot_2, Bween_pool("avg", c1, c2), "avg", c2)
  
  Plot_3 <- Plot_details(pct1, c1, c1$pct.exp, pct1$pct.exp, "% Expressed", "pct", labels[1], 100)

  Plot_4 <- Plot_details(pct2, c2, c2$pct.exp, pct2$pct.exp, "% Expressed", "pct", labels[2], 100)
  Plot_4 <- Plot_Bween(Plot_4, Bween_pool("pct", c1, c2), "pct", c2)
  
  Plot <- Plot_1 + Plot_2 + Plot_3 + Plot_4


  return(Plot)
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

main <- function(){
    
    for (m in 1)
    
    # Create and store a plot for each cluster pool possible clusterpools
    for (subgr_index in subgr) {
      for (id_index in id) {
        # Save the clusterpool that we are working with
        one_cluster_pool <- CPR_new %>%
          filter(SubGroup == subgr_index, id == id_index)
        # Now I should have all the information for one cluster pool
        # creat a temporary list
        tmp_list <- list()
        # Create key names
        key <- paste(subgr_index, id_index, sep = ' ')
        # Next pass the the one_cluster_pool into the Ploting function
        # I might have to relocate.
        # Do the average expression first
        tmp_list[['Average_Expression']] <- Plot_details(one_cluster_pool, 
                                                         cell_roster[[key]], expm1(cell_roster[[key]]$raw_counts), 
                                                         one_cluster_pool$avg.exp, 
                                                         "Average Expression", "avg", key, 12)
        # Then do the percent expressed next
        tmp_list[['Percent_Expressed']] <- Plot_details(one_cluster_pool,
                                                         cell_roster[[key]], expm1(cell_roster[[key]]$raw_counts),
                                                         one_cluster_pool$pct.exp,
                                                         "% Expressed", "pct", key, 100)
        
        # Add the temp list to the plot archive
        plot_arkv[[key]] <- tmp_list
      }
    }
    
    # Make a list of statistics
    ttest_arkv <- list()
    # Now let us build a list of statistics
    for (i in 1:length(cell_roster)) {
      for (j in 1: length(cell_roster)) {
        if (i == j){
          break
        } else if (i != j) {
          
          # creat a temporary list
          tmp_list <- list()
          # Create key names
          key1 <- paste(names(cell_roster)[i], names(cell_roster)[j], sep = 'X')
          key2 <- paste(i, j, sep = '-')
          # Now do average expresson first
          tmp_list[['Average_Expression']] <- Bween_pool("avg", cell_roster[[i]], cell_roster[[j]])
          # Now do percent expressed
          tmp_list[['Percent_Expressed']] <- Bween_pool("pct", cell_roster[[i]], cell_roster[[j]])
          # Might be necessary to switch key 1 with key 2
          tmp_list[['desc_key']] <- key1
          # Add the temp list to the ttest archive
          ttest_arkv[[key2]] <- tmp_list
        }
      }
    }
    # once complete, we must adjust the pvalues
    ttest_arkv <- pvalue_adjust(ttest_arkv) # when I come back I must finish this part
    
    # Now combine graphs for every combination.
    #go through all possible combinations of clusterpools
    for ( i in 1:length(labels(ListByCluster)) ){
      for (j in 1:length(labels(ListByCluster)) ){
        
        if (i == j){
          #break if clusterpools are the same.
          break;
        } else if (i != j){

          #c1 represents first cluster, c2 the second. (there are always only 2 clusters represented in the graph)
          c1 <- ListByCluster[[i]]
          c2 <- ListByCluster[[j]]

          l1 <- labels(ListByCluster)[i]
          l2 <- labels(ListByCluster)[j]

          #call Plot_4_Bar to combine all 4 quadrants with proper design.
          Plot <- Plot_4_Bar(AvgExpPar(c1), PctExpPar(c1), AvgExpPar(c2), PctExpPar(c2), c1, c2, c(l1, l2))
          Plot
          
          #put underlines instead of spaces for the file name
          l1 <- sub(" ", "_", labels(ListByCluster)[i])
          l2 <- sub(" ", "_", labels(ListByCluster)[j])

          #save image onto desktop with custom title.
          save_image(paste('QB_', l1, 'X', l2, sep =""), Plot, width = 5000, height = 5000)
        }
      
      }
    }
    
}

#run this function if you want to load the data
#load_data()


par(mar=c(3,3,1,1))
plot(AllAvg_Exp)


#produce quad barchart with all possible combinations as jpg in a folder in Output.
main()

