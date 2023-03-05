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
