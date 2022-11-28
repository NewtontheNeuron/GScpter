# This script is for removing duplicates or the effect of duplicates

all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_grin.RDS")
e9 <- all_cell_roster[all_cell_roster$cluster == "Excit-09",]
e9[e9$cell.barcode == "SRR6040903_CAGGCCGTGTTN",]

# What levels have duplicates.
# It depends on what groups you have however, for sure you should not have the
# same cell.barcode twice.
# After the filter search for repeats in cell.barcode,
# confirm that it is a duplicate and just remove the duplicate

# Creating a function to get the indicies of the things to be removed
rm_dups <- function(groupeddata) {
  groupeddata[groupeddata$cell.barcode == "SRR6040903_CAGGCCGTGTTN",]
  print(table(groupeddata$cell.barcode))
  groupeddata %>%
    group_by(run) %>%
    count()
}

lbc <- all_cell_roster %>%
  # group by cluster and gene combinations
  group_by(cluster, features.label)# %>%
  #with_groups(NULL, filter)

rm_dups(lbc %>% filter(features.label == "Grin1", cluster == "Excit-09"))
