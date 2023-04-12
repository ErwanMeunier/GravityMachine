library(dplyr)
library(stringr) # truncating strings
library(ggplot2)
library(gridExtra)

setwd("/home/erwan/Documents/GIT/GravityMachine/src/results")

pathBinVar = "resultsBinVar"

# Data must be homogeneous
parsingBinVar <- function(dir){
}

plotBinVar <- function(data){
}

dir = "12 April"
algorithms <- as.list(list.files(paste(pathBinVar,dir,sep="/"), pattern = ".csv", full.names = TRUE, recursive = FALSE))
ref <- read.csv(paste(".",pathBinVar,"ref.csv",sep="/"),sep=";",header=TRUE) # dataframe
#data <- c("ref"=ref, setNames(lapply(algorithms,function(x) read.csv(x,header=TRUE,sep=";")),algorithms))
# # # # # # # # #
refList <- list(ref)
names(refList) = "ref"
# # # # # # # # #
data <- lapply(algorithms,function(x) data.frame(read.csv(x,header=TRUE,sep=";")))
algorithms <- lapply(algorithms, function(x) str_trunc(str_trunc(x,side="left",width=8,ellipsis=""),side="right",width=4,ellipsis=""))
names(data) <- algorithms
# # # # # # # # #
quantities <- colnames(ref)
filenames <- names(data)
result <- do.call(cbind, lapply())