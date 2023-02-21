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
