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
sathy <- read_table(file = "../../../Datasets/Sathyamurthy/GSE103892_Sample_Cell_Cluster_Information.txt") %>%
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