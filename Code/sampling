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

### DOWNLOAD DATA and select desired characteristics ###
#delta
#read tsv file
metadata = as.data.frame(fread('Data/metadata.tsv', drop = c("Type", "Sequence length", 
                                                                   "Host", "Patient age", "Gender", "Clade","Pangolin version", "Variant", 
                                                                   "AA Substitutions", "Submission date", "Is reference?", "Is complete?", 
                                                                   "Is high coverage?", "Is low coverage?", "N-Content", "GC-Content", "Additional location information")))
#filter for variant of concern

metadata_voc <- metadata %>% filter (`Pango lineage` == 'B.1.1.519')

#process metadata

GISAID.processed <- GISAID_process(metadata_voc)

#filter for wanted date 2021-06-01 - 2021-11-30
GISAID.processed <- GISAID.processed %>%
  filter(Date >= as.Date("2020-01-12") & Date <= as.Date("2021-11-30"))
#this creates an object with the epi week
GISAID.processed$epiweek <- as.numeric(get_epiweek(GISAID.processed$Date))

GISAID.processed$epiweek <- GISAID.processed$epiweek + 53

#select only Mexican Sequences

Mexican_sequences <- filter (GISAID.processed, Country == 'Mexico')

#Migration based sub-sampling occurs (excluded due to data permissions)

#load metadata of sequences selected by process

data <- as.data.frame(fread("19nov2021_EPI_ISL_summary_part_1.txt"))

#modify metadata for renaming sequences

data <- select(data, -c(1))

#generation of txt file
df1 <- data%>% unite (hk,'Accession ID': 'Virus name',sep=",")
df2 <- data %>% unite (hk2,'Accession ID': 'Collection date',sep=",")
df3 <- df2 %>% unite (hk3,'hk2': 'Location',sep=",")
df4 <- df3 %>% unite (hk4,'hk3':'Pango lineage',sep=",")
df5 <- bind_cols (df1$hk,df4$hk4)
df6 <- df5 %>% unite (hk,'...1':'...2' ,sep="|")
colnames(df6) <- "EpiID,SeqName,Date,Location,lineage"
df7 <- df6[,1, drop=FALSE]
#to export as text
write.table(df7, file = "19nov2021_EPI_ISL_summary.txt",row.names=FALSE,sep="\t", quote = FALSE)

#Use the csv to copy and paste Accession ID into GISAID to get required sequences
#make sure to point to the correct directory
path <- '/Users/rhysinward/Documents/Mexico_phylogenetics/Github/Data/'
fnames<- dir('/Users/rhysinward/Documents/Mexico_phylogenetics/Github/Data')
fnames<- fnames[grep(".fasta",fnames)]
#remember to rename the file with all the sequences BetaCoV_Wuhan_177seqs_22dec2020 
# create one file with all of the sequences
#need to rename file
combinedFname <- paste("gisaid_hcov-19_2022_01_20_09.fasta",sep="")
pos <- grep(combinedFname,dir('/Users/rhysinward/Documents/Mexico_phylogenetics/Github/Data'),fixed=TRUE)
if (length(pos)==0) {
  for (i in 1:length(fnames)) {
    #seqs <- read.dna( paste(path,fnames[i],sep=""), format="fasta", as.matrix=FALSE)
    #write.dna(seqs,file=paste(path,combinedFname,sep=""),format="fasta",nbcol=-1,colsep="",append=(i>0))
    temp <- readLines(paste(path,fnames[i],sep=""))
    write(temp,file=paste(path,combinedFname,sep=""),append=(i>0))
  }
  print(paste("Written",combinedFname))	
} else {
  print(paste("Already done",combinedFname))
}

#read the sequence information file
info <- read.csv(paste(path,"19nov2021_EPI_ISL_summary.txt",sep=""))
print(info)
#read file and alter the sequence names 
rootname <-	paste(gsub("\\.fas","",combinedFname),"_processedNames",sep="")
pos <- grep(paste(rootname,".fas",sep=""),dir('/Users/rhysinward/Documents/Mexico_phylogenetics/Github/Data'),fixed=TRUE)
if (length(pos)==0) {
  
  seqs <- read.dna(paste(path,combinedFname,sep=""),format="fasta", as.matrix=FALSE)
  taxa <- as.matrix(attributes(seqs)$names)
  isName <- apply(taxa, 1, getEl, ind=3, sep="/")
  epiISL <- apply(taxa, 1, getEl, ind=2, sep="\\|")
  minds  <- match(epiISL, info$EpiID)
  all( epiISL==info$EpiID[minds] )
  dateTxt <- as.matrix(info$Date[minds])
  lineage <- as.matrix(info$lineage[minds])
  decDate <- as.numeric(apply(as.matrix(dateTxt), 1, calcDecimalDate_fromTxt, dayFirst=FALSE, namedMonth=FALSE, sep="-"))
  location<- as.matrix(info$Location[minds])
  country <- apply(as.matrix(location), 1, getEl, ind=1, sep=" / ")
  state   <- apply(as.matrix(location), 1, getEl, ind=2, sep=" / ")
  place   <- apply(as.matrix(location), 1, getEl, ind=3, sep=" / ")
  
  pos     <- which((state=="Kanagawa Prefecture") & is.na(place))
  place[pos] <- "Yokohama"
  print(paste("Changed place of",info$SeqName[pos],"from NA to Yokohama (capital city)"))
  
  newTaxa <- paste(state,isName,epiISL,lineage,dateTxt,sep="|")
  newTaxa <- gsub(" ","_",newTaxa)
  attributes(seqs)$names <- newTaxa
  write.dna(seqs, file=paste(path,rootname,".fas",sep=""), format="fasta", nbcol=-1, colsep="")
  
  newInfo <- cbind(newTaxa,epiISL,isName,dateTxt,country,state,place,decDate)
  colnames(newInfo) <- c("SeqName","EPI_ISL","IsolateName","CollectionDate","Country","State","Place","decDate")
  write.table(newInfo,file=paste(path,rootname,"_infoTbl.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  
} else {
  print("Already done renaming")
}



