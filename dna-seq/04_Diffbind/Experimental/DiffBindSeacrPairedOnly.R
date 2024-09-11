#https://pubmed.ncbi.nlm.nih.gov/26578583/
#https://www.biostars.org/p/9543513/
#https://bioconductor.org/packages/devel/bioc/manuals/DiffBind/man/DiffBind.pdf

setwd("C:\\Users\\hunte\\Desktop\\ChipData")
library(BiocParallel)
register(SerialParam())
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")


#library(ChIPQC)

library(DiffBind)
library(tidyverse)
library(edgeR)
#BiocManager::install("Repitools")
library(Repitools)
#samples <- read.csv('chipSamples2_28_24.csv')

######
#MACS 2 Only Paired End

samplesOnly <- read.csv('chipSamplesOnlyPair.csv')
samplesOnly
dbOnly <- dba(sampleSheet=samplesOnly) #scoreCol=5
dbOnly["Intervals"]
#resultOnly <- dba.count(dbOnly,summits=FALSE)#,summits=20000)
#saveRDS(resultOnly,file="onlyPeakCounts.RDS")
resultOnly=readRDS("onlyPeakCounts.RDS")
resultOnly
testPeaks=dba.peakset(resultOnly)
normOnly <- dba.normalize(resultOnly)
contrastOnly <- dba.contrast(normOnly, minMembers=2)

testPeaks$peaks[3]
testPeaks$samples

#resultOnly <- dba.count(dbOnly,summits = FALSE)#,summits=20000)
#saveRDS(resultOnly,file="onlySummitFalseCounts.RDS")

#https://support.bioconductor.org/p/108094/
#resultOnly$peaks
#normOnly
contrastOnly <- dba.contrast(normOnly,group1=c(1:2,7:8), group2=c(3:4,9:11),
                             name1="K9Me3", name2="K27Ac",minMembers = 2)

contrastOnly <- dba.contrast(normOnly,group1=c(1:2,5:8,12:13), group2=c(3:4,9:11),
                             name1="All_Others", name2="K27Ac")#,minMembers = 2)

contrastOnly <- dba.contrast(normOnly,group1=c(1:2,7:8), group2=c(3:6,9:13),
                             name1="K9Me3", name2="All_Others")#,minMembers = 2)

#saveRDS(contrastOnly,"AcVsAllCHIP.rds")
#saveRDS(contrastOnly,"MeVsAllCHIP.rds")
contrastOnly = readRDS("MeVsAllCHIP.rds")

#contrastOnly <- dba.contrast(normOnly,group1=3:4, group2=9:11,
                             #name1="K9Me3 BF", name2="K9Me3 RVFV",minMembers = 2)

contrastOnly <- dba.contrast(normOnly,group1=3:4, group2=9:11,
                            name1="K27Ac BF", name2="K27Ac RVFV",minMembers = 2)


analizOnly <- dba.analyze(contrastOnly,method=DBA_ALL_METHODS)
par("mar")
par(mar=c(1,1,1,1))
analizOnly
dba.plotVenn(analizOnly,contrast=1,method=DBA_ALL_METHODS)

#analiz <- dba.analyze(contrastOnly,method=DBA_ALL_METHODS)

dba.plotPCA(analizOnly,  attributes=DBA_FACTOR, label=DBA_ID)

plot(analizOnly)

dba.plotVolcano(analizOnly)

dba.plotMA(analizOnly, method=DBA_DESEQ2)

dba.plotMA(analizOnly, bXY=TRUE)

pvals <- dba.plotBox(analizOnly)

res_deseq <- dba.report(analizOnly, method=DBA_EDGER, contrast = 1, th=1)

summary(resultOnly$binding[,3]-resultOnly$binding[,2])

head(n=25,res_deseq)

head(res_deseq,n=9)

library(GenomicRanges)
resDF = annoGR2DF(res_deseq)

sigResDF = head(resDF,n=500,stringsAsFactors=FALSE)


write.csv(sigResDF, "Macs2EdgeRMevsAll.csv", row.names=FALSE)


#saveRDS(sigResDF,"AcVsAllResDF.RDS")

saveRDS(sigResDF,"AcBFVsAcRVFVDF.RDS")

olap.rate <- dba.overlap(dbOnly, mode=DBA_OLAP_RATE)
olap.rate
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')


######
#MACS 2 Strict Paired End


samplesStrict <- read.csv('chipSamplesStrictPair.csv')
dbStrict <- dba(sampleSheet=samplesStrict) #scoreCol=5
#resultStrict <- dba.count(dbStrict)#,summits=20000)
#saveRDS(resultStrict,file="strictPeakCounts.RDS")
resultStrict=readRDS("strictPeakCounts.RDS")
dba.peakset(resultStrict)
normStrict <- dba.normalize(resultStrict)
contrastStrict <- dba.contrast(normStrict, minMembers=2)

#https://support.bioconductor.org/p/108094/
#resultStrict$peaks
#normStrict
contrastStrict <- dba.contrast(normStrict,group1=3:4, group2=9:11,
                             name1="K9Me3 BF", name2="K9Me3 RVFV",minMembers = 2)
analizStrict <- dba.analyze(contrastStrict,method=DBA_ALL_METHODS)
par("mar")
par(mar=c(1,1,1,1))
analizStrict
dba.plotVenn(analizStrict,contrast=1,method=DBA_ALL_METHODS)

analiz <- dba.analyze(contrastStrict,method=DBA_ALL_METHODS)

dba.plotPCA(analizStrict,  attributes=DBA_FACTOR, label=DBA_ID)

plot(analizStrict)

dba.plotVolcano(analizStrict)

dba.plotMA(analizStrict, method=DBA_DESEQ2)

dba.plotMA(analizStrict, bXY=TRUE)

pvals <- dba.plotBox(analizStrict)

res_deseq <- dba.report(analizStrict, method=DBA_EDGER, contrast = 1, th=1)

res_deseq

olap.rate <- dba.overlap(dbStrict, mode=DBA_OLAP_RATE)
olap.rate
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')


######
#SEACR Input Relaxed

samplesInput <- read.csv('seacrNormed.csv')
dbInput <- dba(sampleSheet=samplesInput) #scoreCol=5
resultInput <- dba.count(dbInput,summits=FALSE)#,summits=20000)
#saveRDS(resultInput,file="inputPeakCounts.RDS")
resultInput=readRDS("inputPeakCounts.RDS")
resultInput
dba.peakset(resultInput)
normInput <- dba.normalize(resultInput)
contrastInput <- dba.contrast(normInput, minMembers=2)

#https://support.bioconductor.org/p/108094/
#resultInput$peaks
#normInput
contrastInput <- dba.contrast(normInput,group1=3:4, group2=9:11,
                             name1="K9Me3 BF", name2="K9Me3 RVFV",minMembers = 2)
analizInput <- dba.analyze(contrastInput,method=DBA_ALL_METHODS)
par("mar")
par(mar=c(1,1,1,1))
analizInput
dba.plotVenn(analizInput,contrast=1,method=DBA_EDGER)

dba.plotPCA(analizInput,attributes=DBA_FACTOR, label=DBA_ID)

plot(analizInput)

dba.plotVolcano(analizInput)

dba.plotMA(analizInput, method=DBA_DESEQ2)

dba.plotMA(analizInput, bXY=TRUE)

pvals <- dba.plotBox(analizInput)

res_deseq <- dba.report(analizInput, method=DBA_EDGER, contrast = 1, th=1)

res_deseq

saveRDS(res_deseq,file="SEACRRelaxedPairedOnlyEdgeR_Me_DB.RDS")
res_deseq=readRDS("SEACRRelaxedPairedOnlyEdgeR_Me_DB.RDS")

seqnames = res_deseq[,1][,1]
res_deseq["Conc"]
library(GenomicRanges)
resDF = annoGR2DF(res_deseq)
sigACDF = head(resDF,n=139,stringsAsFactors=FALSE)
#pAdj = res_deseq$FDR

write.csv(sigResDF, "seacrEDGER_Ac.csv", row.names=FALSE)


olap.rate <- dba.overlap(dbInput, mode=DBA_OLAP_RATE)
olap.rate
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')


######
#SEACR No input AUC

# samplesInput <- read.csv('seacrAUC.csv')
# dbInput <- dba(sampleSheet=samplesInput) #scoreCol=5
# resultInput <- dba.count(dbInput)#,summits=20000)
# saveRDS(resultInput,file="aucPeakCounts.RDS")
# resultInput=readRDS("inputPeakCounts.RDS")
# dba.peakset(resultInput)
# normInput <- dba.normalize(resultInput)
# contrastInput <- dba.contrast(normInput, minMembers=2)
# 
# #https://support.bioconductor.org/p/108094/
# #resultInput$peaks
# #normInput
# #testContrast <- dba.contrast(normInput,group1=3:4, group2=12:13,
# #                             name1="K9Me3 BF", name2="K9Me3 RVFV",minMembers = 2)
# analizInput <- dba.analyze(contrastInput,method=DBA_ALL_METHODS)
# par("mar")
# par(mar=c(1,1,1,1))
# analizInput
# dba.plotVenn(analizInput,contrast=1,method=DBA_ALL_METHODS)
# 
# analiz <- dba.analyze(contrastInput,method=DBA_ALL_METHODS)
# 
# dba.plotPCA(analizInput,  attributes=DBA_FACTOR, label=DBA_ID)
# 
# plot(analizInput)
# 
# dba.plotVolcano(analizInput)
# 
# dba.plotMA(analizInput, method=DBA_DESEQ2)
# 
# dba.plotMA(analizInput, bXY=TRUE)
# 
# pvals <- dba.plotBox(analizInput)
# 
# res_deseq <- dba.report(analizInput, method=DBA_DESEQ2, contrast = 1, th=1)
# 
# res_deseq
# 
# olap.rate <- dba.overlap(dbInput, mode=DBA_OLAP_RATE)
# olap.rate
# plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')