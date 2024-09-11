#https://pubmed.ncbi.nlm.nih.gov/26578583/
#https://www.biostars.org/p/9543513/
#https://bioconductor.org/packages/devel/bioc/manuals/DiffBind/man/DiffBind.pdf

setwd("C:\\Users\\hunte\\Desktop\\ChipData")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")


#library(ChIPQC)

library(DiffBind)
library(tidyverse)
library(edgeR)
#samples <- read.csv('chipSamples2_28_24.csv')
samples <- read.csv('chipSamplesCtrled2.csv')
dbObj <- dba(sampleSheet=samples) #scoreCol=5
dbObj
#For less than 3 replicates

#result <- dba(sampleSheet=samples)
#result <- dba.blacklist(dbObj) Need appropriate list....
result <- dba.count(dbObj,summits=20000)
result
dba.peakset(result)
norm_result <- dba.normalize(result)
contrast_norm <- dba.contrast(norm_result)#, minMembers=2)

#https://support.bioconductor.org/p/108094/
dba.contrast()
result$peaks
norm_result
testContrast <- dba.contrast(norm_result,group1=3:4, group2=12:13,
                             name1="K9Me3 BF", name2="K9Me3 RVFV",minMembers = 2)
testAnaliz <- dba.analyze(testContrast,method=DBA_ALL_METHODS)
par("mar")
par(mar=c(1,1,1,1))
testAnaliz
dba.plotVenn(testAnaliz,contrast=1,method=DBA_ALL_METHODS)
res_deseq <- dba.report(testAnaliz, method=DBA_DESEQ2, contrast = 1, th=1)
res_deseq

testContrastAc <- dba.contrast(dbObj,group1=6:7, group2=14:15,
                      name1="Ct1 Cond1", name2="Ct2 Cond2")
testAnalizAc <- dba.analyze(testContrastAc,method=DBA_ALL_METHODS)
dba.plotVenn(testAnalizAc,contrast=1,method=DBA_ALL_METHODS)



testContrastNeg <- dba.contrast(dbObj,group1=6:7, group2=14:15,
                               name1="Ct1 Cond1", name2="Ct2 Cond2")
testAnalizNeg <- dba.analyze(testContrastNeg,method=DBA_ALL_METHODS)
dba.plotVenn(testAnalizNeg,contrast=1,method=DBA_ALL_METHODS)

testContrastIn <- dba.contrast(dbObj,group1=6:7, group2=14:15,
                                name1="Ct1 Cond1", name2="Ct2 Cond2")
testAnalizIn <- dba.analyze(testContrastNeg,method=DBA_ALL_METHODS)
dba.plotVenn(testAnalizNeg,contrast=1,method=DBA_ALL_METHODS)

analiz <- dba.analyze(contrast_norm,method=DBA_ALL_METHODS)

dba.plotPCA(analiz,  attributes=DBA_FACTOR, label=DBA_ID)

plot(analiz)

dba.plotVolcano(analiz)


contrast2 <- dba.contrast(norm_result, categories=DBA_FACTOR, minMembers = 2)
analiz2 <- dba.analyze(contrast2,method=DBA_ALL_METHODS)

dba.plotVenn(analiz,contrast=1,method=DBA_ALL_METHODS)

dba.plotVenn(dbObj, dbObj$masks$MCF7 & dbObj$masks$Responsive)

dba.plotMA(analiz, method=DBA_DESEQ2)

dba.plotMA(analiz, bXY=TRUE)

pvals <- dba.plotBox(analiz)

res_deseq <- dba.report(analiz, method=DBA_DESEQ2, contrast = 1, th=1)

res_deseq
sum(res_deseq$Fold>0)
sum(res_deseq$Fold<0)


out <- as.data.frame(res_deseq)
write.table(out, file="results/Nanog_vs_Pou5f1_deseq2.txt", sep="\t", quote=F, row.names=F)


# Create bed files for each keeping only significant peaks (p < 0.05)

nanog_enrich <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  select(seqnames, start, end)

# Write to file
write.table(nanog_enrich, file="Nanog_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)


pou5f1_enrich <- out %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  select(seqnames, start, end)

# Write to file
write.table(pou5f1_enrich, file="Pou5f1_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)


samples <- read.csv('chipSamplesWInput.csv')
dbInput <- dba(sampleSheet=samples) #scoreCol=5
dbInput
#For less than 3 replicates

#result <- dba(sampleSheet=samples)
#result <- dba.blacklist(dbObj) Need appropriate list....
resultInput <- dba.count(dbInput)
resultInput
normInput <- dba.normalize(resultInput)
contrastInput <- dba.contrast(normInput)#, minMembers=2)
analizInput <- dba.analyze(contrastInput,method=DBA_ALL_METHODS)


dba.plotPCA(analizInput,attributes=DBA_FACTOR, label=DBA_ID)

plot(analizInput)

dba.plotVolcano(analizInput)

reportInput <- dba.report(analizInput)
reportInput

olap.rate <- dba.overlap(dbInput, mode=DBA_OLAP_RATE)
olap.rate
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')

heads(samples)

testContrastCtrl <- dba.contrast(dbInput,group1=3:4, group2=9:10,
                                 name1="BF K27Ac", name2="RVFV K27Ac")


#CTRLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
#LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL




sampleCtrl <- read.csv('chipSamplesCtrledNoInput2.csv')
dbCtrl <- dba(sampleSheet=sampleCtrl) #scoreCol=5
dbCtrl
#For less than 3 replicates

#result <- dba(sampleSheet=samples)
#result <- dba.blacklist(dbObj) Need appropriate list....
resultCtrl <- dba.count(dbCtrl)
resultCtrl
normCtrl <- dba.normalize(resultCtrl)
contrastCtrl <- dba.contrast(normCtrl, minMembers=2)
analizCtrl <- dba.analyze(contrastCtrl,method=DBA_ALL_METHODS)


dba.plotPCA(analizCtrl,attributes=DBA_FACTOR, label=DBA_ID)

plot(analizCtrl)
dba.plotVenn(analizCtrl,contrast=1,method=DBA_ALL_METHODS)


dba.plotVolcano(analizCtrl)

reportCtrl <- dba.report(analizCtrl)
reportCtrl

contrasts <- dba.show(analizCtrl, bContrasts=TRUE)

sampleCtrl

testContrastCtrl <- dba.contrast(dbCtrl,group1=5:6, group2=12:13,
                               name1="BF K27Ac", name2="RVFV K27Ac",contrast=c("Condition","BF","RVFV"))
testContrastCtrl
testAnalizCtrlCon <- dba.analyze(testContrastCtrl,method=DBA_ALL_METHODS,design="~Condition + Factor")
dba.plotVenn(testAnalizCtrlCon,contrast=1,method=DBA_ALL_METHODS)
testAnalizCtrlCon

saveRDS(testAnalizCtrl,"testAnalizCtrl.RDS")
testAnalizCtrl=readRDS("testAnalizCtrl.RDS")
peakReport = dba.report(testAnalizCtrl,th=1,bCounts=TRUE)# with th=1 and bCounts=TRUE
peakReport

dba.plotMA(testAnalizCtrl)

dba.plotPCA(testAnalizCtrl,  attributes=DBA_FACTOR, label=DBA_ID)

plot(testAnalizCtrl)

dba.plotVolcano(testAnalizCtrl)
dba.plotVolcano(analiz)

#dba.plotVenn(dbObj, dbObj$masks$MCF7 & dbObj$masks$Responsive)

dba.plotMA(analiz, method=DBA_DESEQ2)

dba.plotMA(analiz, bXY=TRUE)

pvals <- dba.plotBox(testAnalizCtrl)

res_deseq <- dba.report(analiz, method=DBA_DESEQ2, contrast = 1, th=1)

res_deseq
sum(res_deseq$Fold>0)
sum(res_deseq$Fold<0)






sampleMe <- read.csv('chipSamplesMe.csv')
dbMe <- dba(sampleSheet=sampleMe) #scoreCol=5
dbMe
#For less than 3 replicates

#result <- dba(sampleSheet=samples)
#result <- dba.blacklist(dbObj) Need appropriate list....
resultMe <- dba.count(dbMe)
resultMe
normMe <- dba.normalize(resultMe)
contrastMe <- dba.contrast(normMe,minMembers = 2)#, minMembers=2)
analizMe <- dba.analyze(contrastMe,method=DBA_ALL_METHODS)


dba.plotPCA(analizMe,attributes=DBA_FACTOR, label=DBA_ID)

plot(analizMe)
dba.plotVenn(analizMe,contrast=1,method=DBA_ALL_METHODS)


dba.plotVolcano(analizMe)

reportMe <- dba.report(analizMe)
reportMe

contrasts <- dba.show(analizMe, bContrasts=TRUE)

sampleMe

saveRDS(AnalizMe,"testAnalizMe.RDS")


dba.plotMA(analizMe, method=DBA_DESEQ2)

dba.plotMA(analizMe, bXY=TRUE)

pvals <- dba.plotBox(analizMe)

res_deseq <- dba.report(analizMe, method=DBA_DESEQ2, contrast = 1, th=1)

res_deseq
sum(res_deseq$Fold>0)
sum(res_deseq$Fold<0)


olap.Ctrl.rate <- dba.overlap(dbInput, mode=DBA_OLAP_RATE)
olap.Ctrl.rate
plot(olap.Ctrl.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')
