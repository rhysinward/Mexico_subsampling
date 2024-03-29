GISAID_process <- function (x) {
  #Splitting Location & adding columns 
  by_location <- x %>% separate(Location, into = paste0('Location', 1:4), sep = '/')
  by_location <- select(by_location, c(5:8))
  colnames(by_location) <- c("Continent", "Country","State", "City")
  GISAID.processed <- cbind(x,by_location)
  GISAID.processed$State <- str_trim(GISAID.processed$State, side = c("left"))
  GISAID.processed$State <- str_trim(GISAID.processed$State, side = c("right"))
  GISAID.processed$Country <- str_trim(GISAID.processed$Country, side = c("left")) ##ADDED
  GISAID.processed$Country<- str_trim(GISAID.processed$Country, side = c("right"))
  #Formatting date 
  GISAID.processed$Date.i <- strptime(as.character(GISAID.processed$`Collection date`), "%Y-%m-%d")
  #Remove rows that have incomplete dates e.g.2021-03
  GISAID.processed$Date <- as.character(GISAID.processed$Date.i)
  GISAID.processed <- GISAID.processed[complete.cases(GISAID.processed[ , 28]),]
  #Remove rows that record 'imported' cases
  GISAID.processed <- GISAID.processed[!grepl("Travel history", GISAID.processed$`Additional location information`),]
  #Removing irrelevant columns
  drop <- c("Type", "Collection date", "Location", "Sequence length", 
            "Host", "Patient age", "Gender", "Clade","Pangolin version", "Variant", 
            "AA Substitutions", "Submission date", "Is reference?", "Is complete?", 
            "Is high coverage?", "Is low coverage?", "N-Content", "GC-Content", "Additional location information", "Continent", "City", "Date.i")
  GISAID.processed <- GISAID.processed[,!(names(GISAID.processed) %in% drop)]
}
