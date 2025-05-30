#Reference Publications
#https://pubmed.ncbi.nlm.nih.gov/26578583/
#https://www.biostars.org/p/9543513/
#https://bioconductor.org/packages/devel/bioc/manuals/DiffBind/man/DiffBind.pdf

setwd("C:\\Users\\hunte\\Desktop\\AltChip")
library(BiocParallel)
register(SerialParam())
library(DiffBind)
library(tidyverse)
library(edgeR)
library(GenomicRanges)
library(Repitools)

#BiocManager::install("DiffBind")

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

#Diffbind can use one of two analysis methods for detecting the differential presence of sequences, either DESeq2 or EdgeR
#Both are RNA-seq differential expression libraries and both are reliable, but EdgeR appears to function better with our ChIP-seq software
#EdgeR is the default, but setting 'edger=FALSE' during the function call will set DESeq2 to be used instead

#tableIn allows the user to provide a dataframe object representing the metadata file instead of a file path
#This allows for more flexibility and the use of the metadata sheet for easily selection of samples

runDiffbind = function(sampleFile,samples1,samples2,plotting=FALSE,saving=TRUE,summits=FALSE,sampleName=FALSE,
                       namesOne="Condition1",namesTwo="Condition2",fromFile=FALSE,CSV=TRUE,edger=TRUE,blacklist=c(),tableIn=FALSE){
  if(fromFile){
    contrastOnly=readRDS(sampleFile)
  }
  else{
    if(tableIn){
      samplesOnly = sampleFile
    }
    else{
      samplesOnly <- read.csv(sampleFile)
    }
    #samplesOnly <- read.csv(sampleFile)
    dbOnly <- dba(sampleSheet=samplesOnly) #scoreCol=5
    if(length(blacklist)>0){
      dbOnly=dba.blacklist(dbOnly,blacklist=greylist@regions,greylist=greylist@regions)
    }
    resultOnly <- dba.count(dbOnly,summits=summits)
    if(saving){
      saveRDS(resultOnly,file=paste(namesOne,"vs",namesTwo,"PeakCounts.RDS",sep = ""))
    }
    testPeaks=dba.peakset(resultOnly)
    normOnly <- dba.normalize(resultOnly)
    contrastOnly <- dba.contrast(normOnly,group1=samples1, group2=samples2,
                                 name1=namesOne, name2=namesTwo,minMembers = 2)
  }
  
  analizOnly <- dba.analyze(contrastOnly,method=DBA_ALL_METHODS,bBlacklist = FALSE,bGreylist = FALSE)
  
  
  if(edger){
    res_deseq <- dba.report(analizOnly, method=DBA_EDGER, contrast = 1, th=1)
  }
  else{
    res_deseq <- dba.report(analizOnly, method=DBA_DESEQ2, contrast = 1, th=1)
  }
  
  resDF = annoGR2DF(res_deseq)
  
  if(saving){
    if(CSV){
      write.csv(resDF,paste(namesOne,"vs",namesTwo,".csv",sep=""))#, row.names=FALSE)
    }
    #Change to before to allow different contrasts?
    saveRDS(contrastOnly,paste(namesOne,"vs",namesTwo,".RDS",sep=""))
  }
  
  if(plotting){
    tryCatch(     
      expr = {
        par("mar")
        par(mar=c(1,1,1,1))
        #Error handling?
    
        
        #Good
        dba.plotPCA(analizOnly,  attributes=DBA_FACTOR, label=DBA_ID)
        
        #Good
        plot(analizOnly)
        
        #Error handling?
        
        
        
        
        
        #Error handling?
        if(edger==TRUE){
          dba.plotMA(analizOnly, method=DBA_EDGER)
        }
        else{
          dba.plotMA(analizOnly, method=DBA_DESEQ2)
        }
        dba.plotVenn(analizOnly,contrast=1,method=DBA_ALL_METHODS)
      },
      error = function(e){          
        print("Error: Too few significant peaks too properly plot.")
      }
    )
    
    #?
  }
  return(resDF)
}

#Example Usage


runDiffbind(sampleFile="chipSamplesOnlyPair.csv",samples1 = c(1:3),samples2=c(4:6))

test=readRDS("Condition1vsCondition2.RDS")

newDF=runDiffbind(sampleFile="ChipTestingSamples.csv",samples1 = c(3:4),samples2=c(9:11),namesOne="BF_Ac",namesTwo="RVFV_Ac")

DefaultDBA=readRDS("BF_AcvsRVFV_Ac.RDS")


#Example of selecting samples from larger metadata sheet and using them to run diffbind
allSamps=read.csv("CUT_RUN_Meta_File_MACS2_0.05_keepdup.csv")
day3MeSamps=allSamps[allSamps$Day==3,]
day3MeSamps=day3MeSamps[day3MeSamps$Factor=="H3K9Me3",]


day7MeSamps=allSamps[allSamps$Day==7,]
day7MeSamps=day7MeSamps[day7MeSamps$Factor=="H3K9Me3",]
day7MeSamps=day7MeSamps[day7MeSamps$RiftExperiment==FALSE,]


me_day3=runDiffbind(day3MeSamps,tableIn=TRUE,plotting=TRUE,namesOne="BF_H3K9Me3",namesTwo="SF_H3K9Me3_Day3",samples1=1:3,samples2=4:5)

me_day7=runDiffbind(day7MeSamps,tableIn=TRUE,plotting=TRUE,namesOne="BF_H3K9Me3",namesTwo="SF_H3K9Me3_Day7",samples1=1:3,samples2=4:6)

BFAcSamps=allSamps[allSamps$Factor=="H3K27Ac",]
BFAcSamps=BFAcSamps[BFAcSamps$Condition=="BF",]
BFAcSamps=BFAcSamps[BFAcSamps$Day!=7,]

write.csv(me_day7,"BF_H3K9Me_d7_vs_SF_H3K9Me-DiffBind-out.csv")


BFAc_d1_d3=runDiffbind(BFAcSamps,tableIn=TRUE,plotting=TRUE,namesOne="BF_Ac_D1",namesTwo="BF_Ac_D3",samples1=c(1,3,5),samples2=c(2,4,6),edger = TRUE)

write.csv(BFAc_d1_d3,"BF_Ac_D1vD3DiffbindResults.csv")

BFAcSampsLate=allSamps[allSamps$Factor=="H3K27Ac",]
BFAcSampsLate=BFAcSampsLate[BFAcSampsLate$Condition=="BF",]
BFAcSampsLate=BFAcSampsLate[BFAcSampsLate$Day!=1,]
BFAcSampsLate=BFAcSampsLate[BFAcSampsLate$RiftExperiment=="FALSE",]


BF_Ac_d3_d7=runDiffbind(BFAcSampsLate,tableIn=TRUE,plotting=TRUE,namesOne="BF_Ac_D1",namesTwo="BF_Ac_D3",samples1=c(1:3),samples2=c(4:5),edger = TRUE)
write.csv(BF_Ac_d3_d7,"BF_Ac_D3vD7DiffbindResults.csv")
