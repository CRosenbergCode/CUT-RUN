#Reference Publications
#https://pubmed.ncbi.nlm.nih.gov/26578583/
#https://www.biostars.org/p/9543513/
#https://bioconductor.org/packages/devel/bioc/manuals/DiffBind/man/DiffBind.pdf

#One time Installation
#Run once, then comment out
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")
#BiocManager::install("Repitools")

setwd("C:\\Users\\hunte\\Desktop\\ChipData")
library(BiocParallel)
register(SerialParam())
library(DiffBind)
library(tidyverse)
library(edgeR)
library(GenomicRanges)
library(Repitools)
#samples <- read.csv('chipSamples2_28_24.csv')

runDiffbind = function(sampleFile,samples1,samples2,plotting=FALSE,saving=TRUE,summits=FALSE,sampleName=FALSE,
                       namesOne="Condition1",namesTwo="Condition2",fromFile=FALSE,CSV=TRUE,edger=TRUE){
  #contrastOnly <- dba.contrast(normOnly, minMembers=2)
  
  #resultOnly <- dba.count(dbOnly,summits = FALSE)#,summits=20000)
  #if
  #saveRDS(resultOnly,file="onlySummitFalseCounts.RDS")
  if(fromFile){
    contrastOnly=readRDS(sampleFile)
  }
  else{
    samplesOnly <- read.csv(sampleFile)
    dbOnly <- dba(sampleSheet=samplesOnly) #scoreCol=5
    resultOnly <- dba.count(dbOnly,summits=summits)
    #saveRDS(resultOnly,file="onlyPeakCounts.RDS")
    if(saving){
      saveRDS(resultOnly,file=paste(namesOne,"vs",namesTwo,"PeakCounts.RDS",sep = ""))
    }
    #resultOnly=readRDS("onlyPeakCounts.RDS")
    #resultOnly
    testPeaks=dba.peakset(resultOnly)
    normOnly <- dba.normalize(resultOnly)
    contrastOnly <- dba.contrast(normOnly,group1=samples1, group2=samples2,
                                 name1=namesOne, name2=namesTwo,minMembers = 2)
  }

  analizOnly <- dba.analyze(contrastOnly,method=DBA_ALL_METHODS,bBlacklist = FALSE,bGreylist = FALSE)

  #analiz <- dba.analyze(contrastOnly,method=DBA_ALL_METHODS)

  #pvals <- dba.plotBox(analizOnly)
  if(edger){
    res_deseq <- dba.report(analizOnly, method=DBA_EDGER, contrast = 1, th=1)
  }
  else{
    res_deseq <- dba.report(analizOnly, method=DBA_DESEQ2, contrast = 1, th=1)
  }

  #head(n=25,res_deseq)
  
  #head(res_deseq,n=9)
  
  resDF = annoGR2DF(res_deseq)
  
  #sigResDF = head(resDF,n=25,stringsAsFactors=FALSE)
  if(saving){
    if(CSV){
      write.csv(resDF,paste(namesOne,"vs",namesTwo,".RDS",sep=""), row.names=FALSE)
    }
    saveRDS(resDF,paste(namesOne,"vs",namesTwo,".RDS",sep=""))
  }
  
  if(plotting){
    par("mar")
    par(mar=c(1,1,1,1))
    analizOnly
    dba.plotVenn(analizOnly,contrast=1,method=DBA_ALL_METHODS)
    
    dba.plotPCA(analizOnly,  attributes=DBA_FACTOR, label=DBA_ID)
    
    plot(analizOnly)
    
    dba.plotVolcano(analizOnly)
    
    dba.plotMA(analizOnly, method=DBA_DESEQ2)
    
    dba.plotMA(analizOnly, bXY=TRUE)
    
    plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')
  }
  return(resDF)
}


runDiffbind(sampleFile="chipSamplesOnlyPair.csv",samples1 = c(1:3),samples2=c(4:6))

test=readRDS("Condition1vsCondition2.RDS")

runDiffbind(sampleFile="Condition1vsCondition2PeakCounts.RDS",fromFile = TRUE)

testFile=readRDS("Condition1 vs Condition2")
#
#
#It is highly recommend to save

#runDiffbind = function(sampleCSV="",samples1,samples2,plotting=FALSE,saving=TRUE,summits=FALSE,sampleName=FALSE,
#namesOne="Condition1",namesTwo="Condition2",fromFile=FALSE,CSV=TRUE,edger=TRUE){

#namesOne and namesTwo are not necessary to provide, they are labels
#Using them helps properly keep track of the comparison groups (i.e. RVFV Ac vs. BF Ac)

#'samples1' and 'samples2' are indices which indicate the samples to compare. 
#The numbers are equivalent to rows in the csv or dataframe, excluding headers
#As a reminder, R is 1 indexed so the first object will be 1

#'plotting' simply determines whether or not plots will be displayed
#It is turned off by default

#'saving' allows the user to save intermediate files such as 
#Take care as if the name is the same as a previous file THE OLD FILE WILL BE OVERWRITTEN
#If the file summary file is to be saved, users should do this themselves.

#'fromFile' is a boolean which indicates whether the file provided 
#should be interpreted as a an RDS file holding the counts object
#This is mainly to save time, as the process can
#The samples compared in the contrast can change, but expected peak size
#CANNOT be changed when reading from a file