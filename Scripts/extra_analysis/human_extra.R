# This is a script for further analyses from the human dataset

# Set a key
name <- "Clare"

# Usefull functions
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
p5max <- function(x) {
  round(max(x) * 1.05)
}

# Prepare all_cell_roster if needed
all_cell_roster <- as_tibble(all_cell_roster)

# cluster range plot
Plot <- all_cell_roster %>%
  ggplot(aes(features.label, raw_counts, group = cluster)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.9,
                                              jitter.width = 0,
                                              jitter.height = 0.5),
              aes(color = cluster)) + 
  stat_summary(fun = median, geom = "crossbar",
               position = position_dodge(width = 0.9)) +
  scale_color_viridis_d() +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, p5max(all_cell_roster$raw_counts))) + 
  labs(x = "Gene", y = "Raw log-scaled expression",
       color = "Cluster") +
  ggtitle(wrapper(sprintf("Cluster range plot for human dorsal horn neurons with height jitter- %s", name),
                  width = 45)) +
  #facet_wrap(~subgr) +
  cowplot::theme_cowplot() + 
  theme(legend.position = "bottom",
        axis.title = element_text(size = 15, face = "bold"),
        plot.background = element_rect(fill = "white"))
Plot
save_image("human_nsbgr_hjitter_crplot", Plot, height = 3000, width = 5000)

# table with the number of zeros in every group
znoz <- all_cell_roster %>%
  group_by(features.label) %>%
  summarise(avg.exp = mean(expm1(raw_counts)),
            pct.exp = pct_calc(raw_counts),
            zeros = sum(raw_counts == 0),
            nonzeros = sum(raw_counts != 0),
            variance  = sqrt(sd(raw_counts)))

library(gridExtra)
setwd(sprintf("../Output/%s/", name))
png("DescriptiveStatistics.png")
grid.table(znoz)
dev.off()


# For IASP poster
# Recode the cluster names as ED1 or ID1 etc.
all_cell_roster <- all_cell_roster %>%
  as_tibble() %>%
  mutate(orig.name = cluster,
         cluster = case_when(
           grepl("Ex-", cluster) ~ str_replace(cluster, "Ex-Dorsal-", "ED"),
           grepl("Inh-", cluster) ~ str_replace(cluster, "Inh-Dorsal-", "ID")
         ))
# TODO: restrict or enforce duplicates