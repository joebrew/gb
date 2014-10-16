setwd("~/gb/ebola")

# For first time installation of bioconductor, 
# run the following two lines:
#source("http://bioconductor.org/biocLite.R")
#biocLite()

# For first time installation of GEOquery, run next line:
#biocLite("GEOquery")

# Attach GEOquery library
library(GEOquery)

# download ebola data
#getGEO("GSE24943", destdir="data")

# load into session following download
emat <- getGEO(filename="data/GSE24943_series_matrix.txt.gz")
#softmat <- getGEO(filename="data/GPL9700.soft")

# parse out ID
head(as.character(emat$characteristics_ch1.5), 10)
emat$ID <- sapply(as.character(emat$characteristics_ch1.5), function(x) strsplit(x, ": ")[[1]][2])

# parse out treatment type
emat$tx <- sapply(as.character(emat$characteristics_ch1.3), function(x) strsplit(x, ": ")[[1]][2])

# parse out day of collection
emat$day <- as.numeric(sapply(as.character(emat$characteristics_ch1.4), function(x) strsplit(x, ": ")[[1]][2]))

# have a look
cbind(emat$ID, emat$tx, emat$day)

# assign survival based on table 1 of PMID 21987740
emat$lethal <- 1
emat$lethal[which(emat$ID=="AXX" | emat$ID=="L201-2" | emat$ID=="CH73" | emat$ID=="CL4R")] <- 0

cbind(emat$ID, emat$tx, emat$day, emat$lethal)

# restrict analysis to baseline and 6 days post-infection

# move those with only -8 to baseline (time 0)
table(emat$ID, emat$day)
names(which(table(emat$ID, emat$day)[,2]==0))
emat$day[which(emat$ID %in% names(which(table(emat$ID, emat$day)[,2]==0)) & emat$day==-8)] <- 0
table(emat$ID, emat$day)

# adjust the missing 6-day to be the 3rd day
names(which(table(emat$ID, emat$day)[,3]==0))
emat$day[which(emat$ID %in% names(which(table(emat$ID, emat$day)[,3]==0)) & emat$day==3)] <- 6
table(emat$ID, emat$day)

# subset down
emat <- emat[,which(emat$day %in% c(0,6))]
table(emat$ID, emat$day)

# a look at the data
# observe the missing probes (only 69 are complete, ouch)
temp <- apply(exprs(emat), 1, function(x) sum(is.na(x)))
sum(temp==0)

