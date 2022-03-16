# Start a counter for the total number of comparisons
tcomparisons <- 0
# Function computes a t-test between two cluster pools
# Takes the method = avg or pct, and list by cluster 1 and 2
Bween_pool <- function(method, c1, c2){
  features_no_key <- returnCleanLabelList()
    # Use all_cell_roster and cell_roster inherited from pre_analysis_functions
    GeneStatResults <- data.frame(t = numeric(), df = numeric(), p.value = numeric(), features.label = character())
    GeneStatResults$features.label <- as.character(GeneStatResults$features.label)

    for (i in features_no_key){ 

      if (method == "avg"){
        # The c1 and c2 dataframes need to be:
        # Filtered by the current feature
        # then take the raw counts and expm1 them
        # Then run the ttest
        ClusterStatResults <- t.test(expm1(filter(c1, features.label == i)$raw_counts),
                                     expm1(filter(c2, features.label == i)$raw_counts))
      } else if (method == "pct"){
        ClusterStatResults <- t.test(rnorm(10, 5, 1), rpois(10, 5))
        # Need to fix this later
      }
      # Save the necessary information to Gene stat results
      GeneStatResults <- rbind(GeneStatResults,
                               data.frame(t = ClusterStatResults$statistic,
                                          df = ClusterStatResults$parameter,
                                          p.value = ClusterStatResults$p.value,
                                          features.label = i))
      tcomparisons <<- tcomparisons + 1
    }

    return(GeneStatResults)
}

# Function to apply the pvalue adjustment. Here it is a bonferroni adjustment
# It should take a list of comparison combinations and adjust their pvalues
# based on the tcomparisons counter.
pvalue_adjust <- function(list_of_stats) {
  # Loop through each comparison pair
  for (i in 1:length(list_of_stats)) {
    # Save the current comparison pair to stat
    stat <- list_of_stats[[i]]
    # Loop through each metric/info in the stat object
    for (n in 1:length(stat)) {
      # Exclude the desc_key from getting a p.value not found error
      if (typeof(stat[[n]]) == "list"){
        # Adjust and save the new pvalue
        stat[[n]]$p.value.adj <- stat[[n]]$p.value * tcomparisons
      }
    }
    # Save the new stat object and return
    list_of_stats[[i]] <- stat
  }
  return(list_of_stats)
}

# Function
Plot_Bween <- function(plot, stats, method){
  # Loop throw the row of the ttest dataframe
    for(i in 1:nrow(stats)){
      # Find a value 10% more than the maximum value
      Ext_gene <- filter(plot$data, features.label %in% stats$features.label[i])
      if(method == "avg") {
        value <- max(Ext_gene$avg.exp) + max(Ext_gene$avg.exp) * 0.10
      } else if(method == "pct"){
        value <- max(Ext_gene$pct.exp) + max(Ext_gene$pct.exp) * 0.10
      }
      # Setting the stars based on the adjusted pvalue
      if(stats$p.value.adj[i] < 0.05){
        stars <- "*"
        if(stats$p.value.adj[i] < 0.01){
          stars <- "**"
          if(stats$p.value.adj[i] < 0.001){
            stars <- "***"
          }
        }
        # Change plot if adjusted pvalue is less than 0.05
        plot <- plot + annotate("text", label = stars, x = stats$features.label[i], y = value, size = 6, colour = "black")
      }
    }
    return(plot)
}

#run anova to get p value 
Run_ANOVA <- function (cluster, method, anova.p.val=NULL) {
    if(method == "avg"){
        ANOVA <- aov(expm1(raw_counts) ~ features.label, data = cluster)
    } else if (method == "pct"){
        ANOVA <- aov(expm1(raw_counts) ~ features.label, data = cluster) # This needs to change
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
# Re working plot details to work with CPR
Plot_details <- function (pro_data, clusterpool, clusterpool_exp, method_exp, title, method, label, y_lim){
  features_no_key <- returnCleanLabelList()
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
  
  #add a category title for comparisons between the two clusterpools.
  if (method == "avg"){
    plot <- plot + ggtitle(label) +
      theme(plot.title = element_text(hjust = 0.5, size = 34)) +
      # hard to get standard error of the percent expressed
      geom_errorbar(aes(ymin = avg.lower, ymax = avg.upper),
                    position = position_dodge(width = 0.9), width = 0.50)
  }
  
  return(plot)
  
}

Position_ANOVA <- function(source, method){
  if(method == "avg"){
    pvpos <- max(source$avg.exp) - max(source$avg.exp) * 0.05
  } else if (method == "pct"){
    increment <- 5
    pvpos <- max(source$pct.exp) - max(source$pct.exp) * 0.05
  }
  return(pvpos)
}

mainQBC <- function(CPR){
  id <- returnClusterpool_names()
  subgr <- returnClusterpool_subgroups()
    # Create a list to store all cluster pools and their metrics
    plot_arkv <- list()
    
    # Create and store a plot for each cluster pool possible clusterpools
    for (subgr_index in subgr) {
      for (id_index in id) {
        # Save the clusterpool that we are working with
        one_cluster_pool <- CPR %>%
          filter(subgr == subgr_index, id == id_index)
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
    # go through all possible combinations of clusterpools
    for ( i in 1:length(cell_roster) ){
      for (j in 1:length(cell_roster) ){
        
        if (i == j){
          #break if clusterpools are the same.
          break;
        } else if (i != j){
          
          # Create a comparison key
          comp_key <- paste(i, j, sep = '-')
          # Now at any point in this loop we have
          # The plots for any two cluster pools
          # So all we need to do is put the plots together
          Plot <- plot_arkv[[i]]$Average_Expression +
            Plot_Bween(plot = plot_arkv[[j]]$Average_Expression, stats = ttest_arkv[[comp_key]]$Average_Expression, method = "avg") +
            plot_arkv[[i]]$Percent_Expressed + 
            Plot_Bween(plot = plot_arkv[[j]]$Percent_Expressed, stats = ttest_arkv[[comp_key]]$Percent_Expressed, method = "pct")
          
          #put underlines instead of spaces for the file name
          l1 <- sub(" ", "_", labels(plot_arkv)[i])
          l2 <- sub(" ", "_", labels(plot_arkv)[j])
          #save image onto desktop with custom title.
          save_image(paste('QB_', l1, 'X', l2, sep =""), Plot, width = 5000, height = 5000)
        }
      
      }
    }
    
}