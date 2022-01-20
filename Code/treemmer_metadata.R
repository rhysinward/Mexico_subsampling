# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

#set WD
setwd("/Users/rhysinward/Documents/Mexico_phylogenetics/Github")
# Folder path for results
folres = (path = "./Results/")
# Main functions to run script
files.sources = list.files(path = "./Main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}



#Load packages

library(tidyverse)
library(safejoin)
library(data.table)
library(lubridate)
library(dplyr)
library(zoo) 
library(countrycode) 
library(purrr) 
library(readr) 
library(stringr)   
library(tidyr) 
library(usdata)
library(seqinr)
library(ggplot2)
library(ape)

#Get metadata for treemmer

#names were changed during tree reconstruction so have to extract from tree file
tree <- read.tree("Mex_B_1_1_519_12_07_2021.fasta.treefile")
taxa_spatial_temporal <- data.frame(tree$tip.label)

#spatial_temporal <- read.fasta("updated_519.fas")
#taxa_spatial_temporal <- as.data.frame(as.matrix(attributes(spatial_temporal)$names))
#taxa_split_spatial_temporal <- data.frame(do.call('rbind',strsplit(as.character(taxa_spatial_temporal$V1),'|',fixed = TRUE)))

filter(taxa_spatial_temporal, grepl('Mexico', tree.tip.label))

only_mexico <- filter(taxa_spatial_temporal, grepl('Mexico', tree.tip.label))

rest <- filter(taxa_spatial_temporal, !grepl('Mexico', tree.tip.label))

#newdf1 <- cbind(only_mexico[1], Name=do.call(paste, c(only_mexico, sep="|")))

#newdf2 <- cbind(rest[1], Name=do.call(paste, c(rest, sep="|")))

#newdf1 <- select(newdf1, -c(1))

#newdf2 <- select(newdf2, -c(1))

only_mexico$country <- 'Mexico'
rest$country <- 'Non_mexico'

merged <- rbind (only_mexico,rest)
#merged <- rbind (newdf1,newdf2)

table(merged)

write.table(merged, file = "treemmer_metadata_519.txt",row.names=FALSE,sep=",", quote = FALSE)