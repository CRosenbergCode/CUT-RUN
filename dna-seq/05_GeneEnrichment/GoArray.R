#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#library(readr)
#library("tools")
#library(tidyverse)
setwd("/nfs/home/hogg/GoSplit")
annotFrame=read.table(args[1],sep = "\t",header = TRUE)
goterms=annotFrame$Go.Terms
aeids=annotFrame$Gene.ID

newdf <- data.frame(Gene.ID=character(), 
                    Go.Term=character())

for(i in seq(length(aeids))){
  terms=strsplit(goterms[[i]],split=",")[[1]]
  for(j in terms){
    newdf[nrow(newdf) + 1,] = c(aeids[i],j)
  }
}
write.table(newdf,paste(args[1],".sorted",sep=""),row.names = FALSE,quote=FALSE)


#write.csv(newdf,paste(args[1],".sorted",sep=""))
