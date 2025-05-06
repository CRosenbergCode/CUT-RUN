library(readr)
library(clusterProfiler)
library(tidyverse)
library(fgsea)
library(stringr)
library(tidyverse)
library(topGO)
library(Rgraphviz)
library(GenomicFeatures)
library(ChIPseeker)
library(BiocParallel)

#Necessary for Windows machines
register(SerialParam())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# HOMER peak files should have at minimum 5 columns (separated by TABs, additional columns will be ignored):
#   Column1: Unique Peak ID
# Column2: chromosome
# Column3: starting position
# Column4: ending position
# Column5: Strand (+/- or 0/1, where 0="+", 1="-")
# BED files should have at minimum 6 columns (separated by TABs, additional columns will be ignored)
# Column1: chromosome
# Column2: starting position
# Column3: ending position
# Column4: Unique Peak ID
# Column5: not used
# Column6: Strand (+/- or 0/1, where 0="+", 1="-")


#Takes a diffbind results object and prepares it for usage for homer
DiffbindHomerPrep=function(peakFile,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",promoterOnly=TRUE,enhancerOnly=FALSE,dists=c(-2000,2000),pos=TRUE,neg=FALSE,tableIn=FALSE,returnIDs=TRUE){
  
  tempFile="tempDiffbind.tsv"
  
  testtx=makeTxDbFromGFF(txFile)
  
  if(!tableIn){
    myDiffbind=read.csv(peakFile)
  }
  else{
    myDiffbind=peakFile
  }
  #myDiffbind=read.csv(peakFile)
  slimDiffbind=myDiffbind[c("chr","start","end")]
  
  slimDiffbind["name"]=myDiffbind[[1]]
  slimDiffbind["score"]=0
  slimDiffbind["strand"]="."
  slimDiffbind["signalValue"]=myDiffbind$Fold
  slimDiffbind["pValue"]=-log10(myDiffbind$p.value)
  slimDiffbind["qValue"]=-log10(myDiffbind$FDR)
  slimDiffbind["peak"]=-1
  
  write.table(slimDiffbind,tempFile,sep="\t",row.names=FALSE)
  
  
  myPeaks=annotatePeak(tempFile,tssRegion = dists,TxDb=testtx,verbose=FALSE)
  
  print("PastPeakAnnotation!")
  
  ex_annot = as.data.frame(myPeaks)#myPeaks@anno
  print("The number of rows pre-filtering is")
  print(nrow(ex_annot))
  if(promoterOnly){
    ex_annot=ex_annot[abs(ex_annot$distanceToTSS) <= dists[2],]
  }
  #Best way to incorporate distances?
  #Should this override promoter if TRUE since non-default?
  else if(enhancerOnly){
    ex_annot=ex_annot[abs(ex_annot$distanceToTSS) >= 50000,]
    ex_annot=ex_annot[abs(ex_annot$distanceToTSS) <= 200000,]
    ex_annot=ex_annot[ex_annot$annotation=="Distal Intergenic",]
  }
  if(pos){
    ex_annot=ex_annot[which(ex_annot$signalValue>0),]
  }
  else if (neg){
    ex_annot=ex_annot[which(ex_annot$signalValue<0),]
  }
  print("The number of rows after filtering is")
  print(nrow(ex_annot))


  if(returnIDs){
    unique_annot=ex_annot[!duplicated(ex_annot$geneId),]
    print("After adjusting for uniqueness")
    print(nrow(unique_annot))
    return(unique_annot$geneId)
  }
  else{
    peakformat=ex_annot[c("seqnames","start","end","V4","V5","V6","V7","V8","V9","V10")]
    
    return(peakformat)
  }
}



#Example Usage of function
d1diffbind=read.csv("RVFVd1_H3K27ac_vs_BF-DiffBindout.csv")
d1pos=d1diffbind[which(d1diffbind$Fold>0),]
d1neg=d1diffbind[which(d1diffbind$Fold<0),]

d3diffbind=read.csv("RVFVd3_H3K27ac_vs_BF-DiffBIndout.csv")
d3pos=d3diffbind[which(d3diffbind$Fold>0),]
d3neg=d3diffbind[which(d3diffbind$Fold<0),]


d3diffbind=read.csv("RVFVd3_H3K27ac_vs_BF-DiffBIndout.csv")


sigDay3=d3diffbind[which(d3diffbind$FDR<0.1),]


d3SigPos=DiffbindHomerPrep(sigDay3,pos=TRUE,tableIn=TRUE)
  
write.table(d3SigPos,"D3SigPosPromotersIds.text",row.names = FALSE,col.names = FALSE,quote = FALSE)

d3SigNeg=DiffbindHomerPrep(sigDay3,pos=FALSE,tableIn=TRUE)

write.table(d3SigNeg,"D3SigNegPromotersIds.text",row.names = FALSE,col.names = FALSE,quote = FALSE)

d3Sig=DiffbindHomerPrep(sigDay3,pos=FALSE,neg=FALSE,tableIn=TRUE)


d3Pos=DiffbindHomerPrep("RVFVd3_H3K27ac_vs_BF-DiffBIndout.csv")


write.table(d3Pos,"D3PosPromotersIds.text",row.names = FALSE,col.names = FALSE,quote = FALSE)

d1Neg=DiffbindHomerPrep("RVFVd3_H3K27ac_vs_BF-DiffBIndout.csv",pos=FALSE)

write.table(d1Neg,"D1NegPromotersIds.text",row.names = FALSE,col.names = FALSE,quote = FALSE)


#write.csv(ex_annot,"Ago2D1SFPromoter.csv",row.names = FALSE,quote=FALSE)#,col.names =FALSE)
#write.table(peakformat,"Ago2D1SFPromoter.peakFile",sep="\t",row.names = FALSE,quote=FALSE,col.names =FALSE)
