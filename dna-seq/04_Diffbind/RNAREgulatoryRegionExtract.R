library(DESeq2)
library(ggplot2)
library(pheatmap)
library("vsn")
library(EnhancedVolcano)

setwd("C:\\Users\\hunte\\Desktop\\ChipData")

salmonResults=readRDS("SalmonDESeqResults.RDS")

head(salmonResults)

possibleNegs=head(salmonResults,n=10000)
possibleNegs=tail(possibleNegs,n=9500)
possibleNegs=possibleNegs[abs(possibleNegs$log2FoldChange)<.5,]
nrow(possibleNegs)
possibleNegs=possibleNegs[possibleNegs$lfcSE<.5,]
nrow(possibleNegs)
possibleNegs=possibleNegs[possibleNegs$baseMean>50,]
nrow(possibleNegs)

nrow(possibleNegs)

set.seed(2024)
negRows=sample(seq(nrow(possibleNegs)),1000)
negExamples=possibleNegs[negRows,]

head(negExamples)


posExamples=head(salmonResults,n=1000)
posExamples=posExamples[posExamples$padj<0.1,]
nrow(posExamples)

head(posExamples)

write.csv(posExamples,file="PositiveExamplesRNA7_30.csv")
write.csv(negExamples,file="NegativeExamplesRNA7_30.csv")

library(readr)
wholeBed=read_tsv("AedesGenes.bed",col_names = FALSE)

row.names(posExamples)

wholeBed[11]==row.names(posExamples)

posBed=wholeBed[row.names(posExamples),]

negBed=wholeBed[row.names(negExamples),]


testString=toString(wholeBed[1,10])
substring(testString,4,13)

getGeneName = function(bedFrame){
  nsize=nrow(bedFrame)
  #print(nrow(bedFrame))
  retArr=rep("",nsize)
  for(i in seq(nsize)){
    retArr[i]=substring(toString(bedFrame[i,10]),4,13)
  }
  return(retArr)
}

subsetGeneName = function(bedFrame,namesVec){
  nsize=nrow(bedFrame)
  subsetArr=rep(FALSE,nsize)
  geneNames=getGeneName(bedFrame)
  for(i in seq(nsize)){
    #print(geneNames[i])
    if(toString(geneNames[i]) %in% row.names(namesVec)){
      subsetArr[i]=TRUE
    }
  }
  print(sum(subsetArr))
  return(bedFrame[subsetArr,])
}

testGN=getGeneName(wholeBed)
wholeBed[11]=getGeneName(wholeBed)
head(wholeBed)

#sum(wholeBed[11] %in% row.names(posExamples))
#sum(row.names(posExamples) %in% wholeBed[11])
toString(wholeBed[1,11])
nsize=nrow(wholeBed)
#print(nrow(bedFrame))
retArr=rep(FALSE,nsize)
#for(i in seq(nsize)){
#  if(toString(wholeBed[i,11]) %in% row.names(posExamples)){
#    print("IN")
#  }
  #print("Test")
#}

posBed=subsetGeneName(wholeBed,posExamples)
negBed=subsetGeneName(wholeBed,negExamples)

promStart = lapply(posBed[2], function (x) x-1500)

promEnd = posBed[2]

#head(testCol)

#head(posBed[2])

posBedProm=data.frame(posBed)
posBedProm[2]=promStart
posBedProm[3]=promEnd

promStart = lapply(negBed[2], function (x) x-1500)

promEnd = negBed[2]

negBedProm=data.frame(negBed)
negBedProm[2]=promStart
negBedProm[3]=promEnd
write.table(posBedProm,file="RNAPosProm.txt",quote=FALSE, sep='\t', col.names = FALSE,row.names=FALSE)
write.table(negBedProm,file="RNANegProm.txt",quote=FALSE, sep='\t', col.names = FALSE,row.names=FALSE)