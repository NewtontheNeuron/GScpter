library(Seurat)

# Set wd to the dir of the current script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load('../../Data/human_ariel_data/top_level_new_annotation.rda')
cell_class <- read.csv('../../Data/human_ariel_data/top_level_annotations_metadata.csv')
str(cell_class)
levels(cell_class$x)

# only run the human stuff
RDSfile <- integrated_top_level_obj
clean_neuron_object <- integrated_top_level_obj
rm(integrated_top_level_obj)

# Run both mouse and human
human <- clean_neuron_object
filename <- file.choose()
mouse <- readRDS(filename)
clean_neuron_object <- mouse
gc()

clean_neuron_object <- human
rm(human)

# Find where the clusters are (done outside the loop)
relevant_info <- as.data.frame(clean_neuron_object@meta.data) %>%
  mutate(cells = row.names(as.data.frame(clean_neuron_object@active.ident)),
         id_plus = as.data.frame(clean_neuron_object@active.ident)[[1]])
  
  
clusters_found <- names(which(apply(relevant_info, FUN = function (x) "Excit-01" %in% x, MARGIN = 2) == TRUE))

# this is done inside the loop

CellIndicies <- which(relevant_info[[clusters_found[1]]] %in% c("Excit-01"))


# ---------------------------------
# Looking at sex differences

# Recreate zscoring
zs_calc <- function (x, rm.out = FALSE) {
  if (rm.out){
    ms <- mean(x)
    sds <- sd(x)
    for(i in 1:length(x)) {
      if (x[i] > ms + (4 * sds)){
        x[i] <- NA
        ms <- mean(x, na.rm = TRUE)
        sds <- sd(x, na.rm = TRUE)
      }
    }
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  } else {
    return(x - mean(x)) / sd(x)
  }
}

mad_scaled <- function(x) {
  # First remove those outliers greater than 4* mad
  median_norms <- (x - median(x, na.rm = T)) / mad(x, na.rm = T)
  # Find out which ones are too great and remove them in x
  x[which(median_norms > 3 * mad(x))] <- NA
  # Then recalculate the median norms and return them
  median_norms <- (x - median(x, na.rm = T)) / mad(x, na.rm = T)
  return(median_norms)
}

all_cell_roster <- all_cell_roster %>%
  mutate(ClusterAndSubgroup = paste(id, subgr), 
         age = case_when(
           sample == "Sathyamurthy" ~ "Adult",
           grepl("Rosenberg", sample) ~ "Postnatal",
           grepl("Haring", sample) ~ "Juvenile"),
         clu_sub_age = paste(ClusterAndSubgroup, age))

# lbc
lbc_samp <- all_cell_roster %>%
  # group by cluster and gene combinations
  group_by(cluster, features.label, sample) %>%
  # calculate the average expression and percent expressed
  summarise(avg.exp = mean(expm1(raw_counts)),
            pct.exp = pct_calc(raw_counts),
            ClusterAndSubgroup = paste(id, subgr), 
            age = case_when(
              sample == "Sathyamurthy" ~ "Adult",
              grepl("Rosenberg", sample) ~ "Postnatal",
              grepl("Haring", sample) ~ "Juvenile"),
            id = id,
            clu_sub_age = paste(ClusterAndSubgroup, age)) %>%
  # Ungroup and perform the appropriate scaling
  ungroup(cluster, features.label) %>%
  # This scaling turns to NA any thing 4 standard deviations above the mean
  mutate(avg.exp.scaled = ifelse(avg.exp > mean(avg.exp) + (4 * sd(avg.exp)), NA, log10(avg.exp)))

lbc_samp %>% ggplot(aes(x = cluster, y = log(avg.exp, base = 100),
                        shape = features.label, color = sample)) + geom_point() +
  #geom_line(aes(group = features.label)) +
  theme(axis.text.x = element_text(angle = -45))

lbc_samp %>% ggplot(aes(x = sample, y = log(avg.exp, base = 100),
                        color = cluster)) + geom_point() +
  geom_line(aes(group = cluster)) +
  facet_wrap(~features.label)

lbc_samp %>% ggplot(aes(x = sample, y = avg.exp,
                        color = cluster)) + geom_point() +
  geom_line(aes(group = cluster)) +
  facet_wrap(~features.label)

# CPR
CPR_samp <- all_cell_roster %>%
  # group by clusterpool (id and subgroup) gene combinations
  group_by(id, features.label, sample) %>%
  # perform the following calculations within the groups and store the
  # neccessary information
  summarise(avg.exp = mean(expm1(raw_counts)),
            pct.exp = pct_calc(raw_counts),
            avg.std.err = SE(expm1(raw_counts)),
            avg.lower = avg.exp - avg.std.err,
            avg.upper = avg.exp + avg.std.err) %>%
  # Ungroup the data table and caluclate the appropriate scaling
  ungroup(id, features.label, sample) %>%
  mutate(avg.exp.z.scaled = zs_calc(avg.exp, rm.out = TRUE)) %>%
  mutate(#ClusterAndSubgroup = paste(id, subgr), 
         age = case_when(
           sample == "Sathyamurthy" ~ "Adult",
           grepl("Rosenberg", sample) ~ "Postnatal",
           grepl("Haring", sample) ~ "Juvenile"),
         age_and_id = paste(age, id),
         )#clu_sub_age = paste(ClusterAndSubgroup, age))

# explore
# What I learned the pooling of the pooled dotplot and 
# cluster pooled results relies on the user deciding what 
# they want to pool

# It also relies on how you want things recoded

Plot <- CPR_samp %>%
  filter(age != "Postnatal") %>%
  ggplot(aes(y = features.label, x = age_and_id,
             color = avg.exp.z.scaled, size = pct.exp)) + 
  geom_point() +
  labs(size = '% Expressing', x = "Clusterpools", y = "Gene") +
  scale_size(range = c(0, 20)) +
  scale_color_viridis_c(option = "plasma") + 
  cowplot::theme_cowplot() + 
  theme(axis.title = element_text(size = 20, face = "bold"),
        legend.key.size=unit(1, "line"),
        axis.text.x = element_text(angle = -45, vjust = 0, hjust=0, size=15)) + # changed -45 angle to 0
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15)) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
Plot

Plot <- CPR_samp %>%
  mutate(id = fct_relevel(id, "SDH", "DDH")) %>%
  ggplot(aes(x = age, y = pct.exp, fill = id)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.7) +
  facet_wrap(~features.label, scales = "free_y") +
  scale_fill_hue(h = c(270, 360)) +
  labs(y = "Percent Expressed",
       fill = "Location in dorsal horn",
       x = "Age") +
  #scale_y_continuous(expand = c(0,0)) +
  #cowplot::theme_cowplot() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        legend.key.size = unit(2, "line"),
        legend.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 15),
        panel.background = element_rect(),
        plot.background = element_rect())

Plot

# Reordering
#age = fct_relevel(age,
#                  "Postnatal",
#                  "Juvenile",
#                  "Adult")

# Same thing as above but as a box plot
Plot <- lbc_samp %>%
  mutate(id = fct_relevel(id, "SDH", "DDH")) %>%
  ggplot(aes(x = age, y = avg.exp, fill = id)) +
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.7) +
  facet_wrap(~features.label, scales = "free_y") +
  scale_fill_hue(h = c(270, 360)) + 
  labs(y = "Average Expression",
       fill = "Location in dorsal horn",
       x = "Age") +
  theme(axis.title = element_text(size = 15, face = "bold"),
        legend.key.size = unit(2, "line"),
        legend.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 15),
        panel.background = element_rect(),
        plot.background = element_rect())
Plot

# Now let us look at the distribution differences across the cells
all_cell_roster %>%
  ggplot(aes(x = age, y = expm1(raw_counts), fill = id)) +
  geom_violin(position = position_dodge()) +
  #geom_point() +
  facet_wrap(~features.label, scales = "free_y")

# histograms
Plot <- all_cell_roster %>%
  filter(raw_counts != 0) %>%
  ggplot() +
  geom_histogram(aes(raw_counts, fill = age, color = age), binwidth = 0.1) +
  geom_point(aes(x = raw_counts, y = -40, group = age, color = age), position = position_jitter(height = 10)) +
  coord_cartesian(ylim = c(-50, 450)) +
  facet_wrap(~features.label, scales = "free_y")

# Saving some of my ecplorations
save_image("BoxPlot_Avg_Age", Plot, width = 4700, height = 2300, dpi = 350)

# ---- Explore human data
# Exploring the percent expressed

CPR %>% ggplot(aes(features.label, pct.exp, fill = subgr)) +
  geom_col(position = position_dodge())

Plot <- all_cell_roster %>% ggplot(aes(subgr, raw_counts, fill = sample)) +
  geom_boxplot(position = position_dodge()) +
  facet_wrap(~features.label, scales = "free_y") +
  labs(fill = "Samples",
       y = "Raw log-scaled expression levels",
       x = "Subgroups")

Plot <- CPR %>% ggplot(aes(subgr, log1p(avg.exp), fill = pct.exp)) +
  geom_col(position = position_dodge()) +
  scale_fill_viridis_c() +
  facet_wrap(~features.label, scales = "free_y") +
  labs(fill = "Percent Expressed",
       y = "Log-scaled pooled average expression",
       x = "Subgroups") +
  #scale_y_continuous(expand = c(0,0)) +
  #cowplot::theme_cowplot() +
  theme(axis.title = element_text(size = 10, face = "bold"),
        legend.key.size = unit(1, "line"),
        legend.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 10),
        panel.background = element_rect(),
        plot.background = element_rect())

# Saving some of my ecplorations
save_image("Boxplot_Samples_Human", Plot, width = 4700, height = 2300, dpi = 350)
