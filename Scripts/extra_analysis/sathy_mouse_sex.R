### Sathyamurthy mouse sex and biological sample information

# We will redo the the sex separation and sample separation for the
# sathyamurthy dataset. The first time it was done the new sex labels
# were attached to the metadata and then the metadata was attached
# to all_cell_roster. This time the sex labels themselves will just
# be attached to all_cell_roster.

library(tidyverse)
library(readr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in the table containing the barcodes and their labels
# Also separate the barcodes from the sex/subject label

sathy <- read_delim(file = "../../../Datasets/Sathyamurthy/GSE103892_Sample_Cell_Cluster_Information.txt") %>%
  separate(col = sample_cellbarcode, into = c("sex", "base_barcode"),
           sep = "_")

# Remove the unecessary columns with warnings on import
sathy <- sathy[1:2]
sathy

# Now bring in all_cell_roster
# filter to only sathyamurthy cells
all_cell_roster <- readRDS("../../../Datasets/all_cell_rosters/all_cell_roster_grin.RDS")
all_cell_roster <- as_tibble(all_cell_roster)
all_cell_roster <- filter(.data = all_cell_roster, dataset == "Sathyamurthy")
all_cell_roster

# The barcodes are after the underscore in cell.barcode of all_cell_roster
all_cell_roster <- all_cell_roster %>%
  separate(col = cell.barcode, into = c("cell.barcode.1", "base_barcode"),
           sep = "_")
all_cell_roster

# Check that all the barcodes in all_cell_roster are in sathyamurthy
any(!(unique(all_cell_roster$base_barcode) %in% sathy$base_barcode))
all_cell_roster$base_barcode[which(unique(all_cell_roster$base_barcode) %in% sathy$base_barcode)]

any(sathy$base_barcode %in% all_cell_roster$base_barcode[1])
str_match_all(sathy$base_barcode, all_cell_roster$base_barcode[1])
str_match(sathy$base_barcode, "GTAGCTGG") %>%
  as.vector() %>% .[!is.na(.)]
match_1 <- str_match(sathy$base_barcode, "GTAGCTGG") %>%
  as.vector()
sathy$base_barcode[which(!is.na(match_1))]

# Mark as many as you can
sathy_roster <- all_cell_roster %>%
  left_join(sathy, by = "base_barcode") %>%
  mutate(individual = sex) %>%
  mutate(sex = case_when(
    str_detect(sex, "f|F") ~ "female",
    str_detect(sex, "m|M") ~ "male",
    sex %in% c("rotarod1", "rotarod4", "form3", "form8", "form5", "form10") ~ "female",
    sex %in% c("rotarod2", "rotarod3", "rotarod5",
                     "form1", "form6", "form2", "form7",
                     "form4", "form9") ~ "male"
  ))
sathy_roster
saveRDS(sathy_roster, "../../../Datasets/all_cell_rosters/sathy_grin_acr.RDS")
sathy_roster <- readRDS("../../../Datasets/all_cell_rosters/sathy_grin_acr.RDS")

# Create plots
extra_pool <- list()
extra_pool[["2"]] <- list("id", "sex", "features.label")
extra_pool[["3"]] <- list("sex", "features.label")

# Create the pooled dotplot based on sex and id
# Create cluster pool resutlts
CPR <- createClusterPoolResults(sathy_roster %>%
                                  filter(!is.na(individual)),
                                "2", "z-score")
CPR
CPR <- CPR %>%
  mutate(sex_id = paste(sex, id, sep = " "))

# Then plot the dot plot
Plot <- CPR %>%
  ggplot(aes(y = features.label, x = sex_id,
             color = avg.exp.scaled, size = pct.exp)) + 
  geom_point() +
  labs(size = '% Expressing', x = "Clusterpools", y = "Gene") +
  scale_size(range = c(0, 20)) +
  scale_color_viridis_c(option = "plasma") + 
  cowplot::theme_cowplot() + 
  theme(axis.title = element_text(size = 20, face = "bold"),
        legend.key.size=unit(1, "line"),
        axis.text.x = element_text(angle = -45, vjust = 0, hjust=0, size=15),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
Plot
setwd("../")
args <- c("NMDA_mouse_SDHvDDH", "../../Datasets/neurons_only_2021/clean_neuron_object.RDS")
source("JSON_Handler.r")
save_image("sex_ddh_sdh_adult", Plot, height = 2000, width = 2000)
# Plot only the sex differences for the entire dh
extra_pool[["3"]] <- list("sex", "features.label")
CPR <- createClusterPoolResults(sathy_roster %>%
                                  filter(!is.na(individual)),
                                "3", "z-score")
CPR
Plot <- CPR %>%
  ggplot(aes(y = features.label, x = sex,
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
save_image("sex_dh_adult", Plot, height = 2000, width = 1300)

# Cluster range plot attempt
Plot <- sathy_roster %>%
  filter(!is.na(individual)) %>%
  ggplot(mapping = aes(features.label, raw_counts, group = cluster)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 1,
                                              jitter.width = 0),
              aes(color = cluster)) + 
  stat_summary(fun = median, geom = "crossbar",
               position = position_dodge(width = 1)) +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(expand = c(0,0), limits = c(-1, 5)) + 
  labs(x = "Gene", y = "Raw log-scaled expression",
       color = "Cluster") +
  facet_wrap(~sex + id) +
  cowplot::theme_cowplot() + 
  theme(legend.position = "bottom",
        axis.title = element_text(size = 15, face = "bold"),
        plot.background = element_rect(fill = "white"))
Plot
save_image("sex_adult_crplot", Plot, height = 3000, width = 5000)

# Now look at the dotplot of dh again across individuals
extra_pool[["4"]] <- list("sex", "individual", "features.label")
CPR <- createClusterPoolResults(sathy_roster %>%
                                  filter(!is.na(individual)),
                                "4", "z-score")
CPR
CPR <- CPR %>%
  mutate(concat = paste(sex, individual, sep = " "))

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
        axis.text.x = element_text(angle = -45, vjust = 0, hjust=0, size=15),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
Plot
save_image("individual_dotplot_adult", Plot, height = 3000, width = 5000)

# Each individual is either male or female but we can also look at
# the sdh vs ddh
extra_pool[["5"]] <- list("sex", "individual", "id", "features.label")
CPR <- createClusterPoolResults(sathy_roster %>%
                                  filter(!is.na(individual)),
                                "5", "z-score")
CPR
CPR <- CPR %>%
  mutate(concat = paste(individual, id, sep = " "),
         concat.2 = paste(sex, individual, sep = " "))

Plot <- CPR %>%
  ggplot(aes(y = features.label, x = id,
             color = avg.exp.scaled, size = pct.exp)) + 
  geom_point() +
  labs(size = '% Expressing', x = "Clusterpools", y = "Gene") +
  scale_size(range = c(0, 20)) +
  scale_color_viridis_c(option = "plasma") +
  facet_wrap(~concat.2) +
  cowplot::theme_cowplot() + 
  theme(axis.title = element_text(size = 20, face = "bold"),
        legend.key.size=unit(1, "line"),
        axis.text.x = element_text(angle = -45, vjust = 0, hjust=0, size=15),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
Plot
save_image("invdal_adult_sdh_ddh", Plot, height = 5000, width = 3000)
