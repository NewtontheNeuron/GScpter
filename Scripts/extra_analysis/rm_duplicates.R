# This script is for removing duplicates or the effect of duplicates
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
all_cell_roster <- readRDS("../../../Datasets/all_cell_rosters/all_cell_roster_grin.RDS")
e9 <- all_cell_roster[all_cell_roster$cluster == "Excit-09",]
(e9x <- e9[e9$cell.barcode == "SRR6040903_CAGGCCGTGTTN",])
# There are both SDH and DDH repetitions
# Parts of e9 would be duplicate so are those parts equal? Or which ones are equal?
e9x[1,1:8]
all(e9x[1,1:8] == e9x[8,1:8])
match(e9x[1,1:8] %in% e9x[,1:8])
plyr::match_df(e9x[,1:8], e9x[1,1:8])
e9x$index <- as.numeric(row.names(e9x))
e9x[,-9:-10]
plyr::match_df(e9x[,-9:-10], e9x[1,-9:-10], on = -9)
is.vector(plyr::match_df(e9x[,-9:-10], e9x[1,-9:-10], on = -9)[,9]$index)

# The matches in this case are:
(matchidx <- plyr::match_df(e9x[,-9:-10], e9x[1,-9:-10], on = -9)[,9]$index)
e9x[matchidx]

# for the matches which columns are differnt? And, what are there values?
# We know where we expect differences
e9x[matchidx]$id
# Remove duplicates in the vector
unique(e9x[matchidx]$id)
# I can also expect it to be in subgr because it is another grouping variable.
unique(e9x[matchidx]$subgr)
# Then ask which one is more than one.
if (length(unique(e9x[matchidx]$id)) > 1) {
  idvals <- unique(e9x[matchidx]$id)
}
exists("idvals")
if (length(unique(e9x[matchidx]$subgr)) > 1) {
  subgrvals <- unique(e9x[matchidx]$subgr)
}
exists("subgrvals")
# Then merge duplicates by the ones that have differences
matchidx[2:length(matchidx)]
e9x[-matchidx[2:length(matchidx)]]
# Or you could just change them any way.
e9x$id <- as.list(e9x$id)
e9x[matchidx[1]]$id <- list(idvals)
e9x <- e9x[-matchidx[2:length(matchidx)]]
e9x

# For testing
e9x[nrow(e9x) + 1] <- e9x[7]
e9x[nrow(e9x)]$id <- "VH"
e9x[11]$id <- "DDH"
# subgroup
e9x[nrow(e9x) + 1] <- e9x[nrow(e9x)]
e9x[nrow(e9x)]$subgr <- "Extracitory"

# I have the process now I just need the function
squash_dups <- function(x) {
  x$index <- as.numeric(row.names(x))
  x$id <- as.list(x$id)
  x$subgr <- as.list(x$subgr)
  for (rows in 1:nrow(x)) {
    if (rows > nrow(x)) break
    
    matchidx <- plyr::match_df(x[,-9:-10], x[rows,-9:-10], on = -9)[,9]$index
    if (length(matchidx) == 1) next
    
    if (length(unique(x[matchidx]$id)) > 1) {
      idvals <- unique(x[matchidx]$id)
    }
    if (length(unique(x[matchidx]$subgr)) > 1) {
      subgrvals <- unique(x[matchidx]$subgr)
    }
    
    if(exists("idvals")) {
      x[matchidx[1]]$id <- list(idvals)
    }
    if(exists("subgrvals")) {
      x[matchidx[1]]$subgr <- list(subgrvals)
    }
    
    x <- x[-matchidx[2:length(matchidx)]]
    x$index <- as.numeric(row.names(x))
  }
  x <- x[,1:(ncol(x) -1)]
  return(x)
}
View(squash_dups(all_cell_roster))
# This takes really long so this process done multiple times is not something
# you want to do and it would be better if the vars were stored better in the first
# place