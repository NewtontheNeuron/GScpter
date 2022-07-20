# This is a script for adding the sex label for the sathyamurthy dataset
library(tibble)
library(readr)
# Working with the Sathyamurthy data so let us restrict the meta data
meta_data <- RDSfile@meta.data

meta_data <- meta_data %>%
  rownames_to_column("cells") %>%
  as_tibble() %>%
  filter(dataset == "Sathyamurthy") %>%
  # Separate the cells column into the two components so that
  # we can have the base barcodes for comparison with the
  # sex labeled meta data from sathyamurthy
  separate(col = cells, into = c("run_dbl", "base_barcode"),
           sep = "_")

# Check the origin_ident are the same as the runs dbl
all(meta_data$run_dbl == meta_data$orig.ident)

# import the sex labeled sathyamurthy metadata
sathy <- read_table(file = "C:/Users/no/Documents/Neuroscience - CH BSc/NEUR 4908 - F20 W21/RNA Seq/Data/Sathyamurthy/GSE103892_Sample_Cell_Cluster_Information.txt")
sathy <- read_table(file = "../../../Datasets/sathyamurthy_raw_data/GSE103892_Sample_Cell_Cluster_Information.txt")
# TODO: fix improper import but the first column should be fine
# In fact you can remove the other columns
sathy <- sathy %>%
  separate(col = sample_cellbarcode, into = c("sex_label", "base_barcode"),
           sep = "_")

# What levels are in the sex label column
unique(sathy$sex_label)
# Remove the columns with form or rotarod in the sex_label
sathy <- sathy %>%
  filter(!str_detect(sex_label, "form|rotarod")) %>%
  mutate(sex_label = case_when(
    str_detect(sex_label, "f|F") ~ "female",
    str_detect(sex_label, "m|M") ~ "male",
    sex_label %in% c("rotarod1", "rotarod4", "form3", "form8", "form5", "form10") ~ "female",
    sex_label %in% c("rotarod2", "rotarod3", "rotarod5",
                     "form1", "form6", "form2", "form7",
                     "form4", "form9") ~ "male"
  ))
# What levels are in the sex label column
unique(sathy$sex_label)
# now it is just mail and female

# Remove the other unnecessary columns
sathy <- sathy[1:2]

# There are 6750 cells which is more than the 5510 cells in meta_data

# Now we join the data.frames
# First see if the base_barcodes are present
# TODO: what is going on?
all(meta_data$base_barcode %in% sathy$base_barcode)

meta_sex <- left_join(meta_data, sathy, by = "base_barcode")
any(is.na(meta_sex$sex_label))

# Therefore it could be the removed cells that are not present
# Lets load the data again and keep the cells and allow the join function to
# remove them
sathy <- read_table(file = "C:/Users/no/Documents/Neuroscience - CH BSc/NEUR 4908 - F20 W21/RNA Seq/Data/Sathyamurthy/GSE103892_Sample_Cell_Cluster_Information.txt") %>%
  separate(col = sample_cellbarcode, into = c("sex_label", "base_barcode"),
           sep = "_")

sathy <- sathy[1:2]
sathy

# merge again
any(meta_data$base_barcode %in% sathy$base_barcode)
all(sathy$base_barcode %in% meta_data$base_barcode)

meta_sex <- left_join(meta_data, sathy, by = "base_barcode")
any(is.na(meta_sex$sex_label))

# We still have NAs so this meta data might not be complete.
# There is no way to tell the sex of the form and rotarod models for now
# We can say they are unknown
meta_sex <- meta_sex %>%
  mutate(individual = sex_label) %>%
  mutate(sex_label = case_when(
    #str_detect(sex_label, "form|rotarod") ~ "unknown",
    str_detect(sex_label, "f|F") ~ "female",
    str_detect(sex_label, "m|M") ~ "male",
    sex_label %in% c("rotarod1", "rotarod4", "form3", "form8", "form5", "form10") ~ "female",
    sex_label %in% c("rotarod2", "rotarod3", "rotarod5",
                     "form1", "form6", "form2", "form7",
                     "form4", "form9") ~ "male"
  ))

#### Haring ####

# The Haring barcodes as used by the Ariel's seurat object 
# Uses the a different barcode name from the Haring count matrix
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
  as_tibble(haring_roster)

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
saveRDS(haring_roster, "../../../Datasets/all_cell_rosters/haring_roster_grin.RDS")

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
save_image("sex_ddh_sdh_adult", Plot, height = 2000, width = 2000)
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


# Cluster range plot attempt
haring_roster %>%
  ggplot(mapping = aes(features.label, raw_counts, group = cluster)) + # Human had to use raw_counts instead of
  geom_jitter(position = position_jitterdodge(dodge.width = 1, # expm1(raw counts), which also changed the label
                                              jitter.width = 0), # from unlog-scaled to log-scaled
              aes(color = cluster)) + # Using raw_counts will work with mouse as well and is
  stat_summary(fun = median, geom = "crossbar", # Probably a safe bet, or you could set a parameter in the function
               position = position_dodge(width = 1)) +
  #scale_color_viridis_d(option = "plasma") +
  scale_fill_hue(h = c(180, 270)) +
  scale_y_continuous(expand = c(0,0), limits = c(-1, 12)) + 
  labs(x = "Gene", y = "Raw log-scaled expression",
       color = "Cluster") +
  facet_wrap(~sex) +
  cowplot::theme_cowplot() + 
  theme(legend.position = "bottom",
        axis.title = element_text(size = 15, face = "bold"))
