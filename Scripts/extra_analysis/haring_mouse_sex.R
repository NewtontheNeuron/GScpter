#### Haring ####

# The Haring barcodes as used by the Ariel's seurat object 
# uses the a different barcode name from the Haring count matrix
# There is another Haring csv that can be used to connect the two.

# Load the conversion table.
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
haring_prepath <- "../../../Data/Haring"
#convtable <- read_delim("../../../Datasets/haring_raw_data/GSE103840.txt")
convtable <- read_csv("../../../Datasets/haring_raw_data/sample labels.csv", col_select = 1:2)
convtable

# Load in the count matrix.
library(readxl)
# This step should only be done on a computer that has a lot of ram
countmatrix <- read_excel("../../../Datasets/haring_raw_data/GSE103840_molecule_counts.xlsx", n_max = 2)
countmatrix

# Since age is also a useful value I want to merge the count matrix and convtable
# by cellid and then add the sex and age for each cell into the metadata for the
# haring integrated cells.

# Make the count matrix tidyer
sex_age <- as_tibble(cbind(Title = names(countmatrix), t(countmatrix))) %>%
  slice(-1) %>%
  rename(age = V2, sex = V3) %>%
  mutate(sex = case_when(
    sex == " 1" ~ "F",
    sex == " 0" ~ "M"
  ))
sex_age

# Are all of the cells in the convtable present in the cell id count matrix?
any(!(sex_age$Title %in% convtable$Title))
all(sex_age$Title %in% convtable$Title)

# Merge the two tables
haring <- sex_age %>%
  left_join(convtable, by = "Title")
haring

# Ho many female cells
length(which(haring$sex == "F"))
length(which(haring$sex == "M"))

# Now load in the all_cell_roster
all_cell_roster <- readRDS("../../../Datasets/all_cell_rosters/all_cell_roster_grin.RDS")

haring_roster <- filter(.data = all_cell_roster, dataset == "Haring") %>%
  as_tibble()

# The Accession column in the haring table is similar to the cell.barcode table
# Are all or any of the cells present in both tables?
any(!(unique(haring_roster$cell.barcode) %in% haring$Accession))
all(unique(haring_roster$cell.barcode) %in% haring$Accession)
length(haring$Accession)
length(unique(haring_roster$cell.barcode))
which(!(haring$Accession %in% unique(haring_roster$cell.barcode)))

# Merge the haring table into the haring_roster
names(haring)[4] <- "cell.barcode"
haring_roster <- haring_roster %>%
  left_join(haring, by = "cell.barcode")
haring_roster

# There is now an age x and age y from the two different table coding and conventions
# You can now perform data science on the haring roster data frame and save it
# for later use.
saveRDS(haring_roster, "../../../Datasets/all_cell_rosters/haring_roster_clare.RDS")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../loadLibraries.R")
source("../Pre_analysis_functions.R")
source("../PooledDotPlot.R")
haring_roster <- readRDS("../../../Datasets/all_cell_rosters/haring_roster_grin.RDS")
extra_pool <- list()
extra_pool[["2"]] <- list("id", "sex", "features.label")
extra_pool[["3"]] <- list("sex", "features.label")

# Create the pooled dotplot based on sex and id
# Create cluster pool resutlts
CPR <- createClusterPoolResults(haring_roster, "2", "z-score")
CPR
CPR <- CPR %>%
  mutate(sex = case_when(
    sex == "F" ~ "female",
    sex == "M" ~ "male"
  ),
  sex_id = paste(sex, id, sep = " "))

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
        axis.text.x = element_text(angle = -45, vjust = 0, hjust=0, size=15)) + # changed -45 angle to 0
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=0.5, size=15)) +
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
Plot
setwd("../")
args <- c("NMDA_mouse_SDHvDDH", "../../Datasets/neurons_only_2021/clean_neuron_object.RDS")
source("JSON_Handler.r")
save_image("sex_ddh_sdh_postnatal", Plot, height = 2400, width = 2000)
# No id sep
extra_pool[["3"]] <- list("sex", "features.label")
CPR <- createClusterPoolResults(haring_roster, "3", "z-score")
CPR
CPR <- CPR %>%
  mutate(sex = case_when(
    sex == "F" ~ "female",
    sex == "M" ~ "male"
  ))
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
save_image("sex_dh_postnatal", Plot, height = 2400, width = 1300)

# Cluster range plot attempt
Plot <- haring_roster %>%
  ggplot(mapping = aes(features.label, raw_counts, group = cluster)) + # Human had to use raw_counts instead of
  geom_jitter(position = position_jitterdodge(dodge.width = 1, # expm1(raw counts), which also changed the label
                                              jitter.width = 0), # from unlog-scaled to log-scaled
              aes(color = cluster)) + # Using raw_counts will work with mouse as well and is
  stat_summary(fun = median, geom = "crossbar", # Probably a safe bet, or you could set a parameter in the function
               position = position_dodge(width = 1)) +
  #scale_color_viridis_d(option = "plasma") +
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
save_image("sex_psnl_crplot", Plot, height = 3000, width = 5000)

# Split the haring_roster by individual
haring_roster %>%
  separate(col = Title, into = c("base_barcode", "individual"), sep = "_") %>%
  filter(individual == "A03")
# The ages and sexes are different so that on its own cannot be the
# marker for the individual.
haring_roster %>%
  separate(col = Title, into = c("base_barcode", "individual"), sep = "_") %>%
  filter(individual == "A03", sex == "M", age.y == 19)
# This narrows it down to a few cells and thier genes
