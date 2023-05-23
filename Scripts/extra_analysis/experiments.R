# Experimenting with the rlist package
library(rlist)
library(pipeR)
friends <- list.load("https://renkun-ken.github.io/rlist-tutorial/data/friends.json")

str(friends)
list.search(friends, . == "Ken")

list.search(friends, "Ken" %in% .)
list.search(friends, .[. == "Ken"])
list.search(friends, .[grepl("en", .)])
list.search(friends, .[. == "24"])
list.search(friends, . == "24", classes = "character")
list.search(friends, .[grepl("en", .)], "character", unlist = T)
list.search(friends, .[grepl("en", .)], "character", n = 3, unlist = T)

# List mapping
people <- list.load("https://renkun-ken.github.io/rlist-tutorial/data/sample.json")
str(people)
list.map(people, Age)
list.map(people, "Age")
list.map(people, Name)
people %>>%
  list.filter(Age == 24) %>>%
  list.mapv(.i)
(p <- list.findi(people, Age >= 25, 2))
str(people[[p]])
(p <- list.which(people, "music" %in% Interests))
str(people[p])

# list any
list.any(people, mean(as.numeric(Expertise)) >= 3)
# list class (grouping)
1:10 %>>%
  list.class(. %% 2 == 0)
str(list.class(people, Interests))
str(people)
people %>>%
  list.class(Interests) %>>%
  list.map(. %>>% list.mapv(.name))
people %>>%
  list.class(names(Expertise)) %>>%
  list.map(. %>>% list.mapv(Name))

# Now try it on my own stuff
library(tibble)
library(rlist)
q <- enframe(list(a = 1, b = list(1:4), c = list(4:5), d = list(c(8,2,2,3))))
q
q$value[[1]] <- as.list(q$value[[1]])
list.search(q$value, . == 2)
list.search(q$value, .[. == 2])
list.search(q$value, . == 2, unlist = T)
list.match(q$value, 2 %in% .)


# Trying tha map method from the stack overflow:
# https://stackoverflow.com/questions/42933058/summarizing-with-overlapping-groups-using-dplyr
library(tidyverse)
(test <- data.frame(year = rep(as.character(2014:2016), 2), value = 1:6))
years.groups <- function(start.years, n) {
  ll <- length(start.years)
  start.years <- sort(start.years[-c((ll - (n - 2)):ll)])
  map(as.numeric(as.character(start.years)), function(x) x:(x + (n - 1)))
} # Out puts a list of the things to group

roll.group <- function(vec) {
  vec %>% map_df(~ test %>%
                   filter(year %in% .x) %>%
                   group_by(year = paste(year[which.min(year)], year[which.max(year)], sep = "-")) %>%
                   summarise(value = sum(value)))
}

bind_rows(years.groups(unique(test$year), 2) %>%
            roll.group(),
          years.groups(unique(test$year), 3) %>%
            roll.group(),
          test %>% group_by(year) %>%
            summarise(value = sum(value)))


# Experimenting with lists as factors
mydf<-data.frame(col1=c("a","b"),col2=c("f","j"))
mydf$col1<-as.list(mydf$col1)
mydf$col2<-as.list(mydf$col2)
str(mydf)

# Trying to fix the one column dot plot
# Create a data frame with a cat varible and plot it
library(ggplot2)
library(dplyr)
(td <- data.frame(a = c("a", "b", "b", "a", "b", "c", "c", "a"),
                 b = c(2, 3, 6, 1, 54, 22, 23, 3),
                 c = c(45, 23, 6, 1, 90, 30, 45, 87),
                 d = c("G1", "G1", "G1", "G1", "G2", "G2", "G2", "G2"),
                 e = c("h", "h", "j", "j", "h", "h", "j", "j")))
td$a_e <- paste(td$a, td$e, sep = " ")
td %>%
  ggplot(aes(x = a_e, y = d, color = b, size = c)) +
  geom_point()

# Hmm this works
typeof(td$a)
class(td$a)
str(td$a)
# Add numbers to CPR and plot it.
typeof(CPR$group.label)
class(CPR$group.label)
str(CPR$group.label)

CPR %>%
  ggplot(aes(group.label, features.label,
             color = avg.exp.scaled, size = pct.exp)) +
  geom_point()
# It seems that the factor reordering step is the problem
# I need a more sure fire way of preserving the order in the json
# users can set the order in the json and we will maintain that.

#TODO: after loading ggplot2 and dplyr I get an error at all_cell_roster
# or maybe I loaded the human data with the mouse extrapool


#### dotplot_scales_and_count_normailiation
# Does the log-normalization make the counts more normally distributed?
all_cell_roster %>%
  ggplot() +
  geom_histogram(aes(x = all_cell_roster$raw_counts),
                 fill = "red", binwidth = 0.5) +
  geom_histogram(aes(x = expm1(all_cell_roster$raw_counts)), fill = "blue", alpha = 0.5) +
  scale_y_continuous(limits = c(0, 500))

cpc <- Matrix::colSums(a)
cpg <- Matrix::rowSums(a)

hist(log10(cpc), main = 'counts per cell', col = 'wheat')

cpc <- Matrix::colSums(RDfile@assays$RNA@counts)
cpg <- Matrix::rowSums(RDfile@assays$RNA@counts)

hist(log10(cpc + 1), main = 'counts per cell', col = 'wheat')

# Comparing average expression values
# This portion relies on having RDfile and all_cell_roster
global_size <- 36
seurat_dp <- DotPlot(RDfile, features = c("Grin1", "Grin2a", "Grin2b", "Grin2c",
                                          "Grin2d", "Grin3a", "Grin3b")) +
  coord_flip() +
  #geom_label(aes(label = signif(avg.exp, digits = 2))) +
  labs(y = "Cluster", x = "Gene",
       color = "Avg exp scaled", size = "% Expressing") +
  scale_size(range = c(0, 20)) +
  scale_color_viridis_c(option = "plasma") +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   hjust = 1, size = global_size,
                                   color = "black"),
        axis.text.y = element_text(angle = 0, vjust = 0.5,
                                   hjust = 1, size = global_size,
                                   color = "black"),
        axis.title = element_text(size = global_size, face = "bold"),
        legend.key.size = unit(1.5, "line"),
        legend.text = element_text(size = global_size/1.75),
        legend.title = element_text(size = global_size),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing.x = unit(1.5, "line"),
        plot.background = element_rect(fill = "white"))

#z-score
subdp <- seurat_dp
subdp$data <- subdp$data %>% filter(id %in% unique(all_cell_roster$cluster)) %>%
  mutate(avg.exp.scaled.old = avg.exp.scaled, avg.exp.scaled = zs_calc(avg.exp))
subdp
save_image("SeuDP_GScform_zs", subdp, height = 2600, width = 7600)

#log10
subdp <- seurat_dp
subdp$data <- subdp$data %>% filter(id %in% unique(all_cell_roster$cluster)) %>%
  mutate(avg.exp.scaled.old = avg.exp.scaled, avg.exp.scaled = log10(avg.exp))
subdp
save_image("SeuDP_GScform_log10", subdp, height = 2600, width = 7600)

#log1p
subdp <- seurat_dp
subdp$data <- subdp$data %>% filter(id %in% unique(all_cell_roster$cluster)) %>%
  mutate(avg.exp.scaled.old = avg.exp.scaled, avg.exp.scaled = log1p(avg.exp))
subdp
save_image("SeuDP_GScform_log1p", subdp, height = 2600, width = 7600)

# labels
subdp <- seurat_dp
subdp$data <- subdp$data %>% filter(id %in% unique(all_cell_roster$cluster)) %>%
  mutate(avg.exp.scaled.old = avg.exp.scaled, avg.exp.scaled = log1p(avg.exp))
subdp
save_image("SeuDP_GScform_log1p_lab", subdp, height = 2600, width = 7600)theme_rat <- cowplot::theme_cowplot() +
  theme(
    text = element_text(family = "sans"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.line = element_line(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "transparent"),
    strip.placement = "outside",
    strip.text = element_text(face = "bold"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )





# Human
# This portion relies on having RDfile and all_cell_roster
global_size <- 36
seurat_dp <- DotPlot(RDfile, features = c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2C",
                                          "GRIN2D", "GRIN3A", "GRIN3B"),
                     assay = "RNA",
                     group.by = "new_annotation") +
  coord_flip() +
  geom_label(aes(label = signif(avg.exp, digits = 2))) +
  labs(y = "Cluster", x = "Gene",
       color = "Avg exp scaled", size = "% Expressing") +
  scale_size(range = c(0, 20)) +
  scale_color_viridis_c(option = "plasma") +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   hjust = 1, size = global_size,
                                   color = "black"),
        axis.text.y = element_text(angle = 0, vjust = 0.5,
                                   hjust = 1, size = global_size,
                                   color = "black"),
        axis.title = element_text(size = global_size, face = "bold"),
        legend.key.size = unit(1.5, "line"),
        legend.text = element_text(size = global_size/1.75),
        legend.title = element_text(size = global_size),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing.x = unit(1.5, "line"),
        plot.background = element_rect(fill = "white"))

#z-score
View(seurat_dp$data %>%
  mutate(avg.exp.scaled.new = zs_calc(avg.exp)))
seurat_dp$data %>%
  filter(id %in% unique(all_cell_roster$cluster)) %>%
  mutate(avg.exp.scaled.old = avg.exp.scaled,
         avg.exp.scaled = log10(avg.exp))
new_sdp <- seurat_dp
new_sdp$data <- seurat_dp$data %>%
  filter(id %in% unique(all_cell_roster$cluster)) %>%
  mutate(avg.exp.scaled.old = avg.exp.scaled,
         avg.exp.scaled = log10(avg.exp))
new_sdp

# I think that their stuff is actually a log scale but I will have to confirm
# that by looking into the source code.


save_image("SeuDP_GScform_zs_nokey", seurat_dp, height = 2600, width = 2600)
save_image("SeuDP_GScform_zs", seurat_dp, height = 2600, width = 2600)
save_image("SeuDP_GScform_lab", seurat_dp, height = 2600, width = 2600)

#log10
subdp <- seurat_dp
subdp$data <- subdp$data %>% filter(id %in% unique(all_cell_roster$cluster)) %>%
  mutate(avg.exp.scaled.old = avg.exp.scaled, avg.exp.scaled = log10(avg.exp))
subdp
save_image("SeuDP_GScform_log10", subdp, height = 2600, width = 7600)

#log1p
subdp <- seurat_dp
subdp$data <- subdp$data %>% filter(id %in% unique(all_cell_roster$cluster)) %>%
  mutate(avg.exp.scaled.old = avg.exp.scaled, avg.exp.scaled = log1p(avg.exp))
subdp
save_image("SeuDP_GScform_log1p", subdp, height = 2600, width = 7600)

# labels
subdp <- seurat_dp
subdp$data <- subdp$data %>% filter(id %in% unique(all_cell_roster$cluster)) %>%
  mutate(avg.exp.scaled.old = avg.exp.scaled, avg.exp.scaled = log1p(avg.exp))
subdp
save_image("SeuDP_GScform_log1p_lab", subdp, height = 2600, width = 7600)






# What is equivalent to seurat?
seu_gcp <- lbc %>%
  left_join(b$data %>% mutate(cluster = id, features.label = features.plot),
            by = c("cluster", "features.label"))

all(seu_gcp$avg.exp.x == seu_gcp$avg.exp.y)
all(seu_gcp$pct.exp.x == seu_gcp$pct.exp.y)
# It is really different taking from RNA vs taking from raw
# I might be in some cases understating or overstating the results from my
# visualizations.


# Using the counts per million scale
all_cell_roster$CPM <- run_CPM(all_cell_roster, "CPM")$raw_counts
all_cell_roster$CPT <- run_CPM(all_cell_roster, "countspertotal")$raw_counts
ung_acr <- group_expand(all_cell_roster, overgrouped = "id")
# I had the count * nCount_RNA / 10^6 before
# If I look at the mouse data hardly anything is in the million range
# Why not do counts per thousand?
# It will look the exact same way.
ung_acr %>%
  group_by(id, features.label) %>%
  summarise(pct.exp = pct_calc(raw_counts), avg.exp = mean(CPM)) %>%
  mutate(scale = zs_calc(avg.exp)) %>%
  ungroup(everything()) %>%
  ggplot() +
  geom_point(aes(id, features.label, color = scale, size = pct.exp)) +
  geom_label(aes(id, features.label, label = signif(avg.exp, digits = 2),
                 hjust = 0.5, vjust = -0.5)) +
  scale_color_viridis_c(option = "plasma")


# Can I set a counts per million cutoff for anything too small
# and basically say that it was just undetectable with the selcted assay or 
# conditions.
# Plot it a different way that shows you every once of information
ung_acr %>%
  ggplot() +
  geom_histogram(aes(CPM, fill = id), binwidth = 25) +
  scale_y_continuous(limits = c(0, 10000)) +
  scale_fill_viridis_d() +
  ung_acr %>%
  ggplot() +
  geom_histogram(aes(raw_counts, fill = id), binwidth = 0.1) +
  scale_fill_viridis_d() +
  ung_acr %>%
  ggplot() +
  geom_histogram(aes(nCount_RNA, fill = id), binwidth = 10^4) +
  scale_y_continuous(limits = c(0, 10000)) +
  scale_fill_viridis_d()
# I would love to figure out what that blip is

# Plot how the CPT or CPM changes across 
ung_acr %>%
  ggplot() +
  geom_point(aes(raw_counts, nCount_RNA, size = zs_calc(CPM)),
             color = "grey", alpha = 0.5,
             position = position_jitter(width = 0.15)) +
  ung_acr %>%
  ggplot() +
  geom_point(aes(raw_counts, nCount_RNA, size = zs_calc(CPT)),
             color = "wheat", alpha = 0.5,
             position = position_jitter(width = 0.15))
# Should I unlog the raw counts before calcualting CPM
# CPM is relatively stable
# I wonder what it does to averages
# 


# Staging in ggplot could be valuable to learn
# So lets learn it
mpg
ggplot(mpg, aes(displ, class)) +
  geom_boxplot(outlier.shape = NA)
ggplot(mpg, aes(displ, class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(
    aes(
      label = after_stat(xmax),
      x = stage(displ, after_stat = xmax)
    ),
    stat = "boxplot", hjust = -0.5
  )
ggplot(mpg, aes(displ, class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(
    aes(
      label = after_stat(xmin),
      x = stage(displ, after_stat = xmin)
    ),
    stat = "boxplot", hjust = 1.1
  )
ggplot(mpg, aes(displ, class)) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(
    aes(
      label = after_stat(mean),
      x = stage(displ, after_stat = mean)
    ),
    stat = "boxplot"
  )
ggplot(mpg, aes(class, hwy)) +
  geom_boxplot(aes(colour = class, fill = after_scale(alpha(fill, 0.1))))

ggplot(mpg, aes(class, hwy)) +
  geom_boxplot(aes(fill = stage(hwy, after_scale = alpha(fill, 0.4))))

ggplot(mpg, aes(cty)) +
  geom_density(
    aes(
      group = factor(class),
      color = after_scale(alpha(colour, 0.3)),
      fill = after_stat(sd(x)),
      y = after_stat(count / sum(n[!duplicated(group)]))
    ),
    position = "stack", bw = 1
  ) +
  geom_density(bw = 1)

# Violin plot with the barplot after
index <- 0
pct_rect <- function (x, width, height, location) {
  index <<- index + 1
  dims <- data.frame(mean = mean(x),
                     sd = sd(x),
                     ymin = location - (height / 2),
                     ymax = location + (height / 2)) %>% 
    mutate(seqnum = index,
           xmin = seqnum - (width / 2),
           xmax = seqnum + (width / 2))
  dims
}

ggplot(mpg, aes(class, displ)) +
  geom_violin() +
  stat_summary(
    aes(
      #y = stage(displ, after_stat = 8),
      #label = after_stat(mean),
      fill = after_stat(mean)
    ), geom = "rect", fun.data = ~ pct_rect(.x, width = 0.9,
                                            height = 0.1, location = 8))


get_points <- function(x) {
  subset(x, x < quantile(x, probs = 0.1) | quantile(x, probs = 0.9) < x)
}


ggplot(data = mpg, aes(class, displ)) +
  stat_summary(
    aes(
      fill = after_stat(me)
    ),
    geom="boxplot",
    fun.data = ~ data.frame(ymin = quantile(.x, 0.10),
                            lower = quantile(.x, 0.25),
                            middle = quantile(.x, 0.5),
                            upper = quantile(.x, 0.75),
                            ymax = quantile(.x, 0.90),
                            me = mean(.x),
                            sd = sd(.x))) +
  stat_summary(fun.y = get_points, geom="point")


# C
