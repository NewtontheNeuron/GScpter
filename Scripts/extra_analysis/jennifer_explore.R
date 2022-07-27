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


### Combined male and female dorsal horn comparison between adult and juvenile studies sathyamurthy and haring
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
haring_roster <- readRDS("../../../Datasets/all_cell_rosters/haring_roster_grin.RDS")
sathy_roster <- readRDS("../../../Datasets/all_cell_rosters/sathy_grin_acr.RDS")
extra_pool <- list()
extra_pool[["3"]] <- list("sex", "features.label")

# Recode age as Juvenile or Adult then merge the data tables
haring_roster$age.c <- "Juvenile"
sathy_roster$age.c <- "Adult"
haring_roster$age.c # checking that it added it for all entries

# Full join
new_roster <- full_join(sathy_roster, haring_roster)
new_roster

# There are a lot of empty columns but that is okay they just 
# contain different methods of recording information
# Proceed with Cluster pool results
extra_pool[["5"]] <- list("age.c", "sex", "features.label")
setwd("../")
args <- c("NMDA_mouse_SDHvDDH", "../../Datasets/neurons_only_2021/clean_neuron_object.RDS")
source("loadLibraries.R")
source("Pre_analysis_functions.R")
source("PooledDotPlot.R")
source("JSON_Handler.r")
CPR <- createClusterPoolResults(new_roster %>%
                                  filter(!is.na(sex)), "5", "z-score")
CPR <- CPR %>%
  mutate(sex = case_when(
    sex %in% c("F", "female") ~ "female",
    sex %in% c("M", "male") ~ "male"
  ),
  concat = paste(age.c, sex, sep = " "))
CPR

# Create a dotplot with the concatenated string on the x axis
Plot <- CPR %>%
  ggplot(aes(y = features.label, x = concat,
             color = avg.exp.scaled, size = pct.exp)) + 
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
save_image("sex_age_dh_non_zeros", Plot, height = 2000, width = 2000)

#TODO:
# excluded clusters that have zero means
# What proportion of clusters in a group have a zero mean?
# Where are most of the cluster means?
# Get all the individuals to see how that affects...

# Opening the all_cell_roster
all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_grin.RDS") %>% as_tibble()

z_roster <- all_cell_roster %>%
  filter(dataset == "Zeisel")

z_roster

# Recoding age
all_cell_roster <- all_cell_roster %>%
  mutate(age = case_when(
           dataset == "Sathyamurthy" ~ "Adult",
           grepl("Rosenberg", dataset) ~ "Postnatal",
           grepl("Haring", dataset) ~ "Juvenile"),
         concat.acr = paste(age, subgr))
all_cell_roster
aj_roster <- all_cell_roster %>%
  filter(age != "Postnatal")

# Average expression of each cluster and thier distribution across
# age and laminar location
# Add or remove the raw_counts filter to or not to filter 0 counts
lbc <- all_cell_roster %>%
  filter(!is.na(age)) %>%
  #filter(raw_counts != 0) %>%
  group_by(cluster, features.label, age, id) %>%
  summarise(avg.exp = mean(expm1(raw_counts)))
lbc

Plot <- lbc %>%
  ggplot(aes(avg.exp, features.label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA,
               aes(fill = age),
               position = position_dodge(width = 1)) +
  geom_point(aes(color = cluster, group = age),
              position = position_dodge(width = 1)) +
  labs(x = "Average expression", y = "Gene",
       color = "Cluster", fill = "Age") +
  facet_wrap(~id) +
  scale_x_continuous(expand = c(0, 0)) +
  cowplot::theme_cowplot() + 
  theme(axis.title = element_text(size = 15, face = "bold"),
        plot.background = element_rect(fill = "white"))
Plot
save_image("adult_juv_avgexp_s_ddh_rmz", Plot, height = 1700, width = 4000)
save_image("adult_juv_avgexp_s_ddh", Plot, height = 1700, width = 4000)

# What proportions or each age cluster laminar location group has
# zero raw counts for each of the genes.
zeros <- aj_roster %>%
  group_by(cluster, features.label, age, id) %>%
  summarise(num.z = sum(raw_counts == 0),
            num.nz = sum(raw_counts != 0),
            total = num.z + num.nz,
            pct.nexp = num.z / total * 100,
            pct.exp = num.nz / total * 100)
zeros

Plot <- zeros %>%
  ggplot(aes(pct.exp, features.label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA,
               aes(fill = age),
               position = position_dodge(width = 1)) +
  geom_point(aes(color = cluster, group = age),
              position = position_dodge(width = 1)) +
  labs(x = "Percent Expressed", y = "Gene",
       color = "Cluster", group = "Age") +
  facet_wrap(~id) +
  cowplot::theme_cowplot() + 
  theme(axis.title = element_text(size = 15, face = "bold"),
        plot.background = element_rect(fill = "white"))
save_image("adult_juv_pctexp_s_ddh", Plot, height = 1700, width = 4000)

# Swaping the axes age and id
lbc <- aj_roster %>%
  #filter(!is.na(age)) %>%
  filter(raw_counts != 0) %>%
  group_by(cluster, features.label, age, id) %>%
  summarise(avg.exp = mean(expm1(raw_counts)))
lbc

Plot <- lbc %>%
  ggplot(aes(avg.exp, features.label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA,
               aes(fill = id),
               position = position_dodge(width = 1)) +
  geom_point(aes(color = cluster, group = id),
             position = position_dodge(width = 1)) +
  labs(x = "Average expression", y = "Gene",
       color = "Cluster", fill = "Laminar location") +
  facet_wrap(~age) +
  scale_x_continuous(expand = c(0, 0)) +
  cowplot::theme_cowplot() + 
  theme(axis.title = element_text(size = 15, face = "bold"),
        plot.background = element_rect(fill = "white"))
Plot
save_image("adult_juv_avgexp_s_ddh_rmz", Plot, height = 1700, width = 4000)
save_image("adult_juv_avgexp_s_ddh", Plot, height = 1700, width = 4000)