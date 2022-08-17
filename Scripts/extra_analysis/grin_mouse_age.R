### In this script

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
all_cell_roster <- readRDS("../../Datasets/all_cell_rosters/all_cell_roster_jess.RDS")
all_cell_roster <- all_cell_roster %>%
  mutate(concat_1 = paste(id, subgr), 
         age = case_when(
           age == "Sathyamurthy" ~ "Adult",
           grepl("Rosenberg", age) ~ "Postnatal",
           grepl("Haring", age) ~ "Juvenile"),
         concat_2 = paste(concat_1, age))
all_cell_roster <- as_tibble(all_cell_roster)
aj_roster <- all_cell_roster %>%
  filter(age != "Postnatal")

# for marr and jess
aj_roster <- all_cell_roster

# Create a dot plot with the SDH and DDH separately
# First for the SDH
global_size <- 25 #40
# Function to wrap title text
wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}


extra_pool[["6"]] <- list("features.label", "age", "id")
CPR <- createClusterPoolResults(aj_roster %>%
                                  filter(id == "SDH"), "6", "z-score")
CPR
CPR <- CPR %>%
  mutate(concat = paste(age, id, sep = " "),
         concat = fct_relevel(concat, "Postnatal SDH", "Juvenile SDH",
                              "Adult SDH"))

Plot <- CPR %>%
  ggplot(aes(y = features.label, x = concat,
             color = avg.exp.scaled, size = pct.exp)) + 
  geom_point() +
  labs(size = '% Expressing', x = "Clusterpools", y = "Gene") +
       ggtitle(wrapper("DotPlot - mouse adult SDH and juvenile SDH - Jessica", width = 30)) +
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
save_image("dp_aj_sdh_jess", Plot, height = 2000, width = 1700)

# Now the DDH
CPR <- createClusterPoolResults(aj_roster %>%
                                  filter(id == "DDH"), "6", "z-score")
CPR
CPR <- CPR %>%
  mutate(concat = paste(age, id, sep = " "),
         concat = fct_relevel(concat, "Postnatal DDH", "Juvenile DDH",
                              "Adult DDH"))

Plot <- CPR %>%
  ggplot(aes(y = features.label, x = concat,
             color = avg.exp.scaled, size = pct.exp)) + 
  geom_point() +
  labs(size = '% Expressing', x = "Clusterpools", y = "Gene") +
  ggtitle(wrapper("DotPlot - mouse adult DDH and juvenile DDH - Jessica", width = 30)) +
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
save_image("dp_aj_ddh_jess", Plot, height = 2000, width = 1700)

# SDH vs DDH dotplot
extra_pool[["7"]] <- list("features.label", "id")
CPR <- createClusterPoolResults(aj_roster, "7", "z-score")
CPR

Plot <- CPR %>%
  ggplot(aes(y = features.label, x = id,
             color = avg.exp.scaled, size = pct.exp)) + 
  geom_point() +
  labs(size = '% Expressing', x = "Clusterpools", y = "Gene") +
  ggtitle(wrapper("DotPlot - mouse dorsal horn - Jessica", width = 30)) +
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
save_image("dp_aj_dh_jess", Plot, height = 2000, width = 1300)

# Now see it in the form of a heatmap with both sdh and ddh and
# no zero counts
CPR <- createClusterPoolResults(aj_roster %>%
                                  filter(raw_counts != 0, !(features.label %in% c("Grin2c", "Grin3b"))),
                                "6", "z-score")
CPR
CPR <- CPR %>%
  mutate(concat = paste(age, id, sep = " "),
         concat = fct_relevel(concat, "Postnatal SDH", "Postnatal DDH", "Juvenile SDH", "Juvenile DDH",
                              "Adult SDH", "Adult DDH"))
#CPR$age <- factor(CPR$age, levels =  c("Juvenile", "Adult"))

Plot <- CPR %>%
  ggplot(aes(y = features.label, x = concat,
             fill = avg.exp.scaled)) + 
  geom_tile() +
  labs(x = "Groups", y = "Gene", fill = "Avg exp scaled") +
  #ggtitle(wrapper("Heatmap - mouse adult and juvenile SDH and
  #                DDH the zeros removed- Grin genes", width = 30)) +
  scale_fill_viridis_c(option = "plasma") +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  cowplot::theme_cowplot() + 
  theme(axis.title = element_text(size = global_size, face = "bold"),
        legend.key.size=unit(1, "line"),
        axis.text.x = element_text(angle = -45, vjust = -0.1, hjust = 0.1, size = global_size),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = global_size),
        legend.text = element_text(size = global_size/1.75),
        legend.title = element_text(size = global_size, angle = -90),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))
Plot
save_image("hm_aj_s_ddh_rmz_ntitle_jess", Plot, height = 2000, width = 2000)

# Lets try a heatmap with a counts per million scale
aj_roster
# Create a new roster with the scaled raw counts
aj_roster_cpm <- aj_roster %>%
  mutate(raw_counts = (expm1(raw_counts) / nCount_RNA) * 10^6)
aj_roster_cpm
# Now run the CPR on this with all zeros
aj_roster_cpm %>%
  group_by(features.label, age, id) %>%
  summarise(avg = mean(raw_counts))
CPR <- createClusterPoolResults(aj_roster_cpm, "6", "z-score")
CPR
CPR <- CPR %>%
  mutate(concat = paste(age, id, sep = " "))

Plot <- CPR %>%
  ggplot(aes(y = features.label, x = concat,
             color = avg.exp.scaled, size = pct.exp)) + 
  geom_point() +
  labs(size = '% Expressing', x = "Clusterpools", y = "Gene", color = "Cpm avg exp scaled") +
  ggtitle(wrapper("DotPlot - mouse adult and juvenile SDH and DDH cpmscale- Jessica", width = 30)) +
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
save_image("dp_aj_s_ddh_cpmscale_grin", Plot, height = 2000, width = 1700)

# Lets try removing zeros and putting variance on the size scale
CPR <- createClusterPoolResults(aj_roster %>%
                                  filter(raw_counts != 0), "6", "z-score",
                                sdr = sqrt(sd(expm1(raw_counts))))
CPR
CPR <- CPR %>%
  mutate(concat = paste(age, id, sep = " "),
         zsdr = zs_calc(sdr),
         weightscale = avg.exp * (pct.exp / 100))

Plot <- CPR %>%
  ggplot(aes(y = features.label, x = concat,
             color = weightscale, size = sdr)) + 
  geom_point() +
  labs(size = 'Variance', x = "Groups", y = "Gene", color = "Avg exp z scaled") +
  ggtitle(wrapper("DotPlot - mouse adult and juvenile SDH
                  and DDH w/o zeros but variance on size- Jessica", width = 30)) +
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
save_image("dp_aj_s_ddh_varsize_rmz_grin", Plot, height = 2000, width = 1700)

# Apply the weight scale
CPR <- createClusterPoolResults(aj_roster, "6", "z-score",
                                sdr = sqrt(sd(expm1(raw_counts))))
CPR
CPR <- CPR %>%
  mutate(concat = paste(age, id, sep = " "),
         zsdr = zs_calc(sdr),
         weightscale = avg.exp * (pct.exp / 100),
         zws = zs_calc(weightscale))

Plot <- CPR %>%
  ggplot(aes(y = features.label, x = concat,
             color = zws, size = sdr)) + 
  geom_point() +
  labs(size = 'Variance', x = "Groups", y = "Gene", color = "Avg exp z scaled") +
  ggtitle(wrapper("DotPlot - mouse adult and juvenile SDH
                  and DDH w/o zeros but variance on size- Jessica", width = 30)) +
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
save_image("dp_aj_s_ddh_varsize_rmz_grin", Plot, height = 2000, width = 1700)

# How many non zeros
aj_roster
znoz <- aj_roster %>%
  group_by(features.label, id) %>%
  summarise(avg.exp = mean(expm1(raw_counts)),
            pct.exp = pct_calc(raw_counts),
            zeros = sum(raw_counts == 0),
            nonzeros = sum(raw_counts != 0),
            variance  = sqrt(sd(raw_counts)))

library(gridExtra)
setwd("../Output/Jessica/")
png("cpr_z_nz_var_Aug_15_2022.png", type = "cairo")
grid.table(znoz)
dev.off()
