# This is a script for Jennifer A

all_cell_roster <- as_tibble(all_cell_roster)

sathy_acr <- filter(.data = all_cell_roster, dataset == "Sathyamurthy")

# Change run to cell.barcode before we join
meta_sex <- mutate(.data = meta_sex, cell.barcode = run)
# only merge the columns that you want
meta_sex <- select(.data = meta_sex, last_col(0:2))
# join the dataframe
sathy_acr <- left_join(sathy_acr, meta_sex, by = "cell.barcode")


# Create cluster pool resutlts
CPR_samp <- sathy_acr %>%
  filter(sex_label %in% c("male", "female")) %>%
  # group by clusterpool (id and subgroup) gene combinations
  group_by(id, features.label, sex_label) %>%
  # perform the following calculations within the groups and store the
  # neccessary information
  summarise(avg.exp = mean(expm1(raw_counts)),
            pct.exp = pct_calc(raw_counts),
            avg.std.err = SE(expm1(raw_counts)),
            avg.lower = avg.exp - avg.std.err,
            avg.upper = avg.exp + avg.std.err) %>%
  # Ungroup the data table and caluclate the appropriate scaling
  ungroup(everything()) %>%
  mutate(avg.exp.z.scaled = zs_calc(avg.exp),
         sex_id = paste(sex_label, id))


# Then plot the dot plot
Plot <- CPR_samp %>%
  ggplot(aes(y = features.label, x = sex_id,
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

save_image("sex_ddh_sdh_adult", Plot, height = 2000, width = 2000)

# It should be possible to do a cluster range plot here
sathy_acr %>%
  filter(sex_label %in% c("male", "female")) %>%
  ggplot(aes(features.label, raw_counts, group = cluster)) +
  geom_point(position = position_dodge(width = 1),
             size = 2.5, mapping = aes(color = cluster)) +
  stat_summary(fun = mean, geom = "crossbar", position = position_dodge(width = 1)) +
  facet_wrap(~sex_label) +
  scale_color_viridis_d() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Gene", y = "Raw log-scaled expression",
       color = "Clusters") +
  theme(title = element_text(size = 20),
        legend.position = "none",
        axis.text.x = element_text(angle = 0, color = "black", size = 20),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 20, face = "bold"),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

length(which(sathy_acr$sex_label == "male"))
length(which(sathy_acr$sex_label == "female"))

#### between animal ####
# Create cluster pool resutlts
CPR_samp <- sathy_acr %>%
  # group by clusterpool (id and subgroup) gene combinations
  group_by(id, features.label, individual) %>%
  # perform the following calculations within the groups and store the
  # neccessary information
  summarise(avg.exp = mean(expm1(raw_counts)),
            pct.exp = pct_calc(raw_counts),
            avg.std.err = SE(expm1(raw_counts)),
            avg.lower = avg.exp - avg.std.err,
            avg.upper = avg.exp + avg.std.err) %>%
  # Ungroup the data table and caluclate the appropriate scaling
  ungroup(everything()) %>%
  mutate(avg.exp.z.scaled = zs_calc(avg.exp),
         sex_id = paste(individual, id))

# Then plot the dot plot
Plot <- CPR_samp %>%
  ggplot(aes(y = features.label, x = sex_id,
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


#### finding umap ####
CleanNeuro_UMAP <- RDSfile@reductions$CleanNeuro_UMAP@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var ="cell.barcode") %>%
  as_tibble() %>%
  right_join(all_cell_roster, by = "cell.barcode")

CleanNeuro_UMAP %>%
  ggplot(aes(CleanNeuro_1, CleanNeuro_2, color = raw_counts)) +
  geom_point() +
  #scale_color_viridis_c(option = "plasma") +
  facet_wrap(~features.label)

# It would be good to do a separate pca and a separate umap
