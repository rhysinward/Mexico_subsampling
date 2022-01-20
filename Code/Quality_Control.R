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


#select for only sequences which have been deemed 'good' by NextClade 
#load results of QC

nextclade <- as.data.frame(fread("nextclade/output/nextclade.tsv"))

#select for sequences that weren't deemed good 

aligned_seq <- read.fasta("nextclade/data/sars-cov-2/sequences.aln.fasta")

taxa <- as.data.frame(as.matrix(attributes(aligned_seq)$names))

nextclade_remove <- filter(nextclade, qc.overallStatus != 'good')

species.to.remove <- nextclade_remove$seqName

vec.names<-unlist(lapply(strsplit(names(aligned_seq), ";"), function(x)x[length(x)]))

vec.tokeep <-which(! vec.names %in%  species.to.remove)

length(vec.tokeep)

write.fasta(sequences=aligned_seq[vec.tokeep], names=names(aligned_seq)[vec.tokeep], file.out="updated_519.fas")

