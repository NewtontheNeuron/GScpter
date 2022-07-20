# This is a script to be translated to an R markdown
# that will show an analysis of the mouse and human datasets
# related to the Phox2a gene

# Hypothesis
# https://www.nature.com/articles/s41598-021-97105-w
# Alsulaiman et al. (2021) suggest that Phox2a is preferentially expressed in 
# lamina 1 projection neurons of the spinal cord. These projection neurons
# have synapses on the ipsilateral side but send signals to the brain on the contralateral side
# In the mouse and human single cell/nuclear RNA seq datasets, we have neurons 
# from the dorsal horn of the spinal cord that could be projection of interneurons.
# Therefore, our research questions are could we use this Phox2a gene to identify lamina
# 1 projection neurons within the clusters, how much Phox2a is expressed in the mouse
# SDH vs DDH, are there more excitatory or inhibitory projection neurons, are there
# any genes that are co-expressed with Phox2a - can we use it to estimate which neurons
# in general are located in lamina 1, what is the distribution of the Phox2a gene across
# clusters and cells, how many cells in the SDH and DDH highly express Phox2a, can we re-cluster the cells based on Phox2a expression.

# Based on the information and research questions I hypothesize that:
# - Phox2a expression patterns will have a range/variance that suggests the presence of
# a different subtype of dorsal horn neurons
# - There will be higher average Phox2a expression in the mouse SDH than the mouse DDH.
# - It is possible that some clusters may have more Phox2a expression than others because
# they were formed by looking for differentially expressed genes.
# - Phox2a is more highly expressed in mouse and human dorsal horn cells labeled as 
# excitatory rather than those labeled as inhibitory.
# - We will need to look more into the protein product of Phox2a and its roles
# know if it can be used to preferentially identify other lamina 1 neurons.
# - Phox2a more higly expressed than NMDAR subunit genes


# The analysis

# We will begin by using the pipeline on 4 different analysis projects:
# mouse SDH vs DDH excitatory vs inhibitiory

# Based on the dot plot and the cluster range plot it seems that not a lot
# of cells express the Phox2a gene. Could it be more localized in the cytoplasm.
# In fact it seems that the DDH excitatory neurons have higher percent expressed.
# Three clusters by eye have higher percent expressed than any of the other clusters.
# Those clusters are Excit 21, 23 and 29. According to Russ et al. (2021) the highly
# expressed genes in these clusters are: Lmx1b, Zfhx3, Nms, Lypd1 for Excit 21; Lmx1b,
# Nfib, Cep112, Cdh23, and Satb1 for Excit 23; and Onecut2 and Pmfbp1 for Excit 29.
# We will have to do some extra research to see if these genes relate to Phox2a.

# If we look at Seqseek and look for other Phox genes as well what patterns emerge.
# Niether Phox2a nor Phox2b are expressed highly in a large amount of cells.

# mouse SDH vs DDH

# Based on the above results we can skip this step, and consider if Phox2a is present in human
# lamina 1 projection neurons and not mice lamina 1 projection neurons.

# human excitatory vs inhibitiory

# A quick visit to the shiny app from Yadav et al. (2022) shows that there
# there are very little transcript levels of Phox2a/b in the any of the cells.

# human all pooled together

# We can skip this step


# Discussion/conclusion
# Unfortuantley, non of the hypotheses are true/valid. Both the mouse and human data sets show that Phox2a/Phox2b transcript levels is very low across the majority of neuronal clusters in the SDH or DDH. Very low transcript levels means that it would be hard for us to draw conclusions because a lot of factors can contribute to cells having low transcript levels. Could it be that Phox2a is only expressed when a projection signal is being sent, or is it not localized in the nucleus? Low transcript levels could be a result of strong regulation or the technique itself. For instance, what was the Cre or BAC protocol doing that distinctly differes from the single cell RNA seq experiments. Albeit we are comparing two methods and there could be flaws with both. There are a few mouse DDH clusters that have a decent percent expressed of Phox2a, which could be investigated further in the future. The only conclusion that can be drawn is that it is unlikely that Phox2a can be used with these single cell/nuclear RNA seq mouse and human datasets to identify projection neurons.


#---
# Update on June 3, 2022
# Get all_cell_roster from the pipline
all_cell_roster <- as_tibble(all_cell_roster)
all_cell_roster

# How many cells in total
total <- all_cell_roster %>%
  filter(features.label == "Phox2a") %>%
  group_by(id, subgr, cluster) %>%
  count(name = "totalcells")

# I want to know how many cells express Phox2a and I also need to know
# where they are.
all_cell_roster %>%
  filter(features.label == "Phox2a", raw_counts > 0) %>%
  group_by(id, subgr, cluster) %>%
  count(name = "numberofcells") %>%
  left_join(total, by = c("id", "subgr", "cluster"))
