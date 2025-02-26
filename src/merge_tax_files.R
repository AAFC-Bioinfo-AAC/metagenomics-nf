#!/usr/bin/env Rscript
library("tidyr")
setwd("./")
# Read all species level files into lists

data_list <- lapply(list.files(path="Species/",
                               pattern="*species.summary.tsv",
                               full.names = T),
                    FUN=function(x) {read.delim(x, header=T)} )

# Merge list into long form dataframe

data_long <- Reduce(function(x,y) merge(x,y,all=TRUE), data_list)

# Some  cleanup
data_long$reads <- NULL
data_long$file <- gsub(".out", "", data_long$file)

# Convert long form to wide dataframe
data_wide <- spread( data_long, file, percent)

# Order by taxon_name
data_wide <- data_wide[ order(data_wide$taxon_name), ]

#Write results to csv
write.csv( data_wide, file="kaiju_merged_species.csv", row.names = F)
