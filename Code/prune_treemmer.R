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

#keep taxa selected for by treemmer

data <- as.data.frame(fread("Mex_B_1_1_519_12_07_2021.fasta.treefile_trimmed_list_X_4000",sep = "",header = FALSE))

setwd("/Users/rhysinward/Documents/Mexico_phylogenetics")

aligned_QC_seq <- read.fasta("updated_519.fas")

QC_taxa <- as.data.frame(as.matrix(attributes(aligned_QC_seq)$names))

#taxa names were changed when tree was made so need to marry it up 

#extract ID for sequences

sequence_ID <-  data.frame(do.call('rbind',strsplit(as.character(QC_taxa$V1),'|',fixed = TRUE)))

#extract ID for tree

tree_ID <-  data.frame(do.call('rbind',strsplit(as.character(data$V1),'_',fixed = TRUE)))

tree_ID <-  data.frame(do.call('rbind',strsplit(as.character(tree_ID$X3),'|',fixed = TRUE)))

species.to.keep <- tree_ID$X1

#vec.names<-unlist(lapply(strsplit(names(aligned_QC_seq), ";"), function(x)x[length(x)]))

vec.to.keep <-which(sequence_ID$X3 %in%  species.to.keep)

length(vec.to.keep)

write.fasta(sequences=aligned_QC_seq[vec.to.keep], names=names(aligned_QC_seq)[vec.to.keep], file.out="updated_519_treemmer.fas")

#remove outliars from tempest 
#make file with epiID of sequences wanting to be removed
info <- read.csv(paste("sequence_to_remove_tempest_519.csv",sep=""))
print(info)
combinedFname <- paste("updated_519_treemmer.fas",sep="")
seqs <- read.dna(paste(combinedFname,sep=""),format="fasta", as.matrix=FALSE)
taxa <- as.data.frame(as.matrix(attributes(seqs)$names))
isName <- apply(taxa, 1, getEl, ind=1, sep="\\|")
epiISL <- apply(taxa, 1, getEl, ind=2, sep="\\|")
taxa$minds  <- as.numeric(match(isName, info$EpiID))
dim(taxa)
taxa$minds <- as.numeric(taxa$minds)
taxa <- filter(taxa, minds > 0)
#load package
library(seqinr)
#load file containing sequences
data<-read.fasta("updated_519_treemmer.fas")
#create a vector containing species names: these are the last part of the string
vec.names<-unlist(lapply(strsplit(names(data), ";"), function(x)x[length(x)]))
#find names to keep: indices which are not in the species to remove
species.to.remove<-taxa$V1
vec.tokeep<-which(! vec.names %in%  species.to.remove)
length(vec.tokeep)
#write the final output
write.fasta(sequences=data[vec.tokeep], names=names(data)[vec.tokeep], file.out="updated_519_treemmer_tempest.fas")
