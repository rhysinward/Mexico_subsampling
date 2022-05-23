#Merida workshop 

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

#set WD
setwd("~/rhys/Mexico")
# Folder path for results
folres = (path = "./Results/")
# Main functions to run script
files.sources = list.files(path = "./Main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./Main/", files.sources[i]), collapse = ''))
}

#Load packages

library(tidyverse)
library(safejoin)
library(data.table)
library(lubridate)
library(dplyr)
library(stringr)   
library(tidyr) 
library(seqinr)
library(ggplot2)
library(ape)
library(purrr)

### PROCESS GISAID METADATA ###

#Load metadata

metadata = as.data.frame(fread('14Dic21_Mex_DeltaQual.tsv',drop = c("Type", "Sequence length", 
                                                                                "Host", "Patient age", "Gender", "Clade","Pangolin version", "Variant", 
                                                                                "AA Substitutions", "Submission date", "Is reference?", "Is complete?", 
                                                                                "Is high coverage?", "Is low coverage?", "N-Content", "GC-Content", "Additional location information",
                                                                                "Sampling strategy","Patient status","Last vaccinated","Passage","Specimen","Additional host information",
                                                                                "Lineage")))

#Process metadata
GISAID.processed <- dplyr :: select(metadata, c(1,2,3,4,5))

### Add epiweek ###

GISAID_final <- GISAID.processed %>%
  mutate(Date = as.Date(`Collection date`, "%d/%m/%Y")) %>%
  filter(Date >= "2021-04-21", Date <= "2021-11-30") %>%
  mutate(epiweek = epiweek(Date))

#Plot raw data

raw_data_plot <- GISAID_final %>%
  group_by(epiweek) %>%
  dplyr :: summarize(count = n())
colnames(raw_data_plot) <- c("Epiweek", "count")

ggplot(raw_data_plot, aes(x= Epiweek, y = count)) + geom_bar(stat = 'identity', position = 'stack') +
  labs(x="Epiweek",y="Total Number of Sequences") + theme_bw() 

#Select top 6 lineages only 

length(table(GISAID_final$`Pango lineage`))

top_lineges <- GISAID_final %>%
  group_by(`Pango lineage`) %>%
  dplyr :: summarize(count = n())
colnames(raw_data_plot) <- c("lineage", "count")

top_6_lineages <- top_lineges[top_lineges$count >= top_lineges$count[order(top_lineges$count, decreasing=TRUE)][6] , ]

ggplot(top_6_lineages, aes(x= `Pango lineage`, y = count)) + geom_bar(stat = 'identity', position = 'stack') +
  labs(x="Epiweek",y="Total Number of Sequences") + theme_bw() 

gisaid_lineage_sampled <- filter(GISAID_final, `Pango lineage` %in% top_6_lineages$`Pango lineage`)

### Proportional subsampling ###
#Will make proportional to the number of cases

#get cases from - https://datos.covid-19.conacyt.mx/#DownZCSV

cases <- read.csv("Casos_Diarios_Estado_Nacional_Confirmados_20220515.csv", header = TRUE)

cases <- filter(cases, nombre == 'Nacional')

# first remember the names
n <- cases$nombre
# transpose all but the first column (name)
df.aree <- as.data.frame(t(cases[,c(-1:-3)]))
colnames(df.aree) <- n
df.aree$date <- factor(row.names(df.aree))

#remove row.names

rownames(df.aree) <- NULL

#remove X in name of date

df.aree$date <- gsub("X", "", df.aree$date)

#Calculate epiweek 

table(GISAID_final$Date)

#Remove incorect dates as first genome was found april 21st
#https://tecreview.tec.mx/2021/07/07/en/covid-19-variants-of-concern-circulating-in-mexico/

final_cases <- df.aree %>%
  mutate(Date = as.Date(date, "%d.%m.%Y")) %>%
filter(Date >= "2021-04-21", Date <= "2021-11-30") %>%
         mutate(epiweek = epiweek(Date)) %>%
  group_by(epiweek) %>%
  dplyr :: summarize(count = sum(Nacional))

#plot cases

ggplot(final_cases, aes(x= epiweek, y = count)) + geom_bar(stat = 'identity', position = 'stack') +
  labs(x="Epiweek",y="Total Number of Cases") + theme_bw() 

#Make proportional to 1 

final_cases$proportion <- final_cases$count/(sum(final_cases$count))

#want 2838 sequences

final_cases$ideal_proportion <- final_cases$proportion*2838

#get sequences

#remove week 17 as no sequences

final_cases <- filter(final_cases, epiweek != 17)

set_seed = 7

proportional_sequences <- gisaid_lineage_sampled %>% 
  group_split(epiweek) %>% 
  map2_dfr(c(final_cases$ideal_proportion), ~ slice_sample(.x, n = .y))

#plot proportional data

proportional_data_plot <- proportional_sequences %>%
  group_by(epiweek) %>%
  dplyr :: summarize(count = n())
colnames(proportional_data_plot) <- c("Epiweek", "count")

ggplot(proportional_data_plot, aes(x= Epiweek, y = count)) + geom_bar(stat = 'identity', position = 'stack') +
  labs(x="Epiweek",y="Total Number of Sequences") + theme_bw() 

proportional_data_plot <- proportional_sequences %>%
  group_by(`Pango lineage`) %>%
  dplyr :: summarize(count = n())
colnames(proportional_data_plot) <- c("Lineage", "count")

ggplot(proportional_data_plot, aes(x= Lineage, y = count)) + geom_bar(stat = 'identity', position = 'stack') +
  labs(x="Lineage",y="Total Number of Sequences") + theme_bw() 

#want to match the Accession ID to original metadata 

final_proportional_sample_orginal_metadata <- metadata[metadata$`Accession ID` %in% proportional_sequences$`Accession ID`, ]

#Match to genomic data

seq <- read.fasta("DELTA_BLANCA_aln.fasta")

taxa <- as.data.frame(as.matrix(attributes(seq)$names))

taxa <- as.data.frame(taxa[!duplicated(taxa$V1), ])

remove <- filter(metadata, ! `Accession ID` %in% (final_proportional_sample_orginal_metadata$`Accession ID`))

species.to.remove <- remove$`Virus name`

vec.names<-unlist(lapply(strsplit(names(seq), ";"), function(x)x[length(x)]))

vec.names <- data.frame(do.call('rbind',strsplit(as.character(vec.names),'|',fixed = TRUE)))

vec.tokeep <-which(! vec.names$X1 %in%  species.to.remove)

length(vec.tokeep)

write.fasta(sequences=seq[vec.tokeep], names=names(seq)[vec.tokeep], file.out="Data/proportion_data.fasta")

#This is used for renaming of genomic data

#write table
write.table(final_proportional_sample_orginal_metadata, file = "Data/11may2022_EPI_ISL_summary_part_1.txt",row.names=TRUE,sep=",", quote = TRUE)
#want csv so can look up Accession ID to get required sequences
write.csv(final_proportional_sample_orginal_metadata, file = "Data/11may2022_EPI_ISL_summary_part_1.csv",row.names=TRUE, quote = TRUE)

### RENAME SEQUENCES ###

#load in data

data <- as.data.frame(fread("Data/11may2022_EPI_ISL_summary_part_1.txt"))

#remove unwanted columns

data <- dplyr ::select(data, -c(1,7,8,9))

#generation of txt file
df1 <- data%>% unite (hk,'Accession ID': 'Virus name',sep=",")
df2 <- data %>% unite (hk2,'Collection date',sep=",")
df3 <- df2 %>% unite (hk3,'hk2': 'Location',sep=",")
df4 <- df3 %>% unite (hk4,'hk3':'Pango lineage',sep=",")
df5 <- bind_cols (df1$hk,df4$hk4)
df6 <- df5 %>% unite (hk,'...1':'...2' ,sep=",")
colnames(df6) <- "EpiID,SeqName,Date,Location,lineage"

#to export as text
write.table(df6, file = "Data/11may2022_EPI_ISL_summary.txt",row.names=FALSE,sep="\t", quote = FALSE)

#Use the csv to copy and paste Accession ID into GISAID to get required sequences

#make sure to point to the correct directory
path <- "~/rhys/Chile/Data/"
fnames <- dir(path)
fnames <- fnames[grep(".fasta",fnames)]
# create one file with all of the sequences
#need to rename file
combinedFname <- paste("proportion_data.fasta",sep="")
pos <- grep(combinedFname,dir("~/rhys/Chile/Data/"),fixed=TRUE)

#read the sequence information file
info <- read.csv(paste(path,"11may2022_EPI_ISL_summary.txt",sep=""))
info$delta <- 'Delta_sub'
print(info)
#read file and alter the sequence names 
rootname <-	paste(gsub("\\.fas","",combinedFname),"_beastNames",sep="")
pos <- grep(paste(rootname,".fas",sep=""),dir(),fixed=TRUE)
if (length(pos)==0) {
  
  seqs <- read.dna(paste(path,combinedFname,sep=""),format="fasta", as.matrix=FALSE)
  taxa <- as.matrix(attributes(seqs)$names)
  isName <- apply(taxa, 1, getEl, ind=3, sep="/")
  epiISL <- apply(taxa, 1, getEl, ind=1, sep="\\|")
  minds  <- match(epiISL, info$SeqName)
  all( epiISL==info$epiISL[minds] )
  dateTxt <- as.Date(info$Date,"%d/%m/%Y")[minds]
  location<- as.matrix(info$Location[minds])
  epiISL <- as.matrix(info$EpiID[minds])
  country   <- apply(taxa, 1, getEl, ind=2, sep="/")
  lineage <- as.matrix(info$lineage[minds])
  delta <- as.matrix(info$delta)
  
  newTaxa <- paste(country,location,isName,epiISL,lineage,delta,dateTxt,sep="|")
  newTaxa <- gsub(" ","_",newTaxa)
  attributes(seqs)$names <- newTaxa
  write.dna(seqs, file=paste(path,rootname,".fas",sep=""), format="fasta", nbcol=-1, colsep="")
  
  newInfo <- cbind(newTaxa,epiISL,isName,dateTxt,country)
  colnames(newInfo) <- c("SeqName","EPI_ISL","IsolateName","CollectionDate","Country","decDate")
  write.table(newInfo,file=paste(path,rootname,"_infoTbl.txt",sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  
} else {
  print("Already done renaming")
}

#run through QC pipeline found here -https://clades.nextstrain.org/

#load results of QC

nextclade <- as.data.frame(fread("nextclade.tsv"))

#select for sequences that weren't deemed good 

aligned_seq <- read.fasta("sequences.aln.fasta")

taxa <- as.data.frame(as.matrix(attributes(aligned_seq)$names))

taxa <- as.data.frame(taxa[!duplicated(taxa$V1), ])

nextclade_remove <- filter(nextclade, qc.overallStatus != 'good')
species.to.remove <- nextclade_remove$seqName

vec.names<-unlist(lapply(strsplit(names(aligned_seq), ";"), function(x)x[length(x)]))

vec.tokeep <-which(! vec.names %in%  species.to.remove)

length(vec.tokeep)

write.fasta(sequences=aligned_seq[vec.tokeep], names=names(aligned_seq)[vec.tokeep], file.out="Data/final_proportional_data_qc.fas")

#load fasta

final_seqs <- read.fasta("Data/final_proportional_data_qc.fas")
final_seqs_name <- as.data.frame(as.matrix(attributes(final_seqs)$names))
final_seqs_name <- data.frame(do.call('rbind',strsplit(as.character(final_seqs_name$V1),'|',fixed = TRUE)))

#plot 

proportional_data_plot <- final_seqs_name %>%
  group_by(X6) %>%
  dplyr :: summarize(count = n())
colnames(proportional_data_plot) <- c("Date", "count")

proportional_data_plot <- proportional_data_plot %>% group_by(week = cut(as.Date(Date), "week", start.on.monday = FALSE)) %>% 
  summarise(count = sum(count)) 


ggplot(proportional_data_plot, aes(x= as.Date(week), y = count)) + geom_bar(stat = 'identity', position = 'stack') +
  labs(x="Date",y="Total Number of Sequences") + theme_bw()  +
  scale_x_date(date_breaks = "1 month", date_labels = "%Y %b %d")

proportional_data_plot <- final_seqs_name %>%
  group_by(X5) %>%
  dplyr :: summarize(count = n())
colnames(proportional_data_plot) <- c("Lineage", "count") 

ggplot(proportional_data_plot, aes(x= Lineage, y = count)) + geom_bar(stat = 'identity', position = 'stack') +
  labs(x="Lineage",y="Total Number of Sequences") + theme_bw()



