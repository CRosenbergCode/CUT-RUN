setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#Required libraries
library(ChIPseeker)
library(clusterProfiler)
library(AnnotationDbi)
library(GenomicFeatures)

#Inputs
#peakFileList: This is the 
#Optional: sampNames is a vector which is empty by default. If the vector is not empty and is equal to the number of samples, the names will be applied to samples for both printing and plotting
#Optional: The path to an annotation file, generally a GFF. This should contain different parts of genes (UTRs, exons, etc) in addition to the whole gene body
#Optional: If verbose is set to true, a summary of each peaks annotation will be printed.

#Output
#An array of annotated peak files, 1 per input sample

plotPeakDistances = function(peakFileList,sampNames=c(),verbose=TRUE,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff"){
  samplefiles <- as.list(samplefiles)
  
  testtx=makeTxDbFromGFF(txFile,format="gff3")
  
  
  peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=testtx, 
                         tssRegion=c(-1000, 1000), verbose=FALSE)
  
  
  if(verbose){
    print(peakAnnoList)
  }
  if(length(names)>0){
    names(peakAnnoList)=sampNames
  }
  par("mar")
  par(mar=c(1,1,1,1))
  plotAnnoBar(peakAnnoList)
  plotDistToTSS(peakAnnoList, title="Distribution of Peaks\n relative to TSS")
  return(peakAnnoList)
}

#Example usage
samplefiles=c("C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewSEACRPeaks\\BF_D7_Ac_Rep1.relaxed.bed","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewSEACRPeaks\\BF_D7_Ac_Rep2.relaxed.bed","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewSEACRPeaks\\RVFV_D7_Ac_Rep1_1.relaxed.bed","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewSEACRPeaks\\RVFV_D7_Ac_Rep1_2.relaxed.bed","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewSEACRPeaks\\RVFV_D7_Ac_Rep2.relaxed.bed")


exampleDistances=plotPeakDistances(samplefiles)
plotAnnoBar(exampleDistances)
plotDistToTSS(exampleDistances, title="Distribution of Peaks\n relative to TSS")


#Required libraries
library(ggplot2)
library(idr)
library(readr)
library(GenomicRanges)


#IDR can only be calculated by comparison of two files at a time

#Inputs
#peakFile1 is 
#peakFile2 is a second file with the same format as one
#Optional: Column with ranking metric to use (7, the (), by default)
#

#Outputs

getIDR = function(peakFile1,peakFile2,){
  df1 <- read_delim(peakFile1, delim="\t", col_names=FALSE)
  df2 <- read_delim(peakFile2, delim="\t", col_names=FALSE)
  
  peak1 <- GRanges(df1$X1, IRanges(df1$X2, df1$X3), score=df1$X8)
  peak2 <- GRanges(df2$X1, IRanges(df2$X2, df2$X3), score=df2$X8)
  length(peak1)
  length(peak2)
  head(peak1)
  fo <- findOverlaps(peak1, peak2)
  length(peak1)
  
  length(peak2)
  
  length(fo)
  
  table(duplicated(from(fo)))
  
  table(duplicated(to(fo)))
  
  fo <- as.data.frame(fo)
  fo <- fo[!duplicated(fo$queryHits) & !duplicated(fo$subjectHits),]
  
  y1 <- peak1$score[fo[,1]]
  y2 <- peak2$score[fo[,2]]
  plot(y1, y2, cex=.1)
  
  plot(log10(y1), log10(y2), cex=.1)
  
  plot(rank(-y1), rank(-y2), cex=.1)
  
  dat <- cbind(log10(y1), log10(y2))
  #dat <- dat[sample(nrow(dat),2000),]
  system.time({ 
    res <- est.IDR(dat, mu=3, sigma=1, rho=.9, p=.5)
  })
  
  df <- data.frame(rep1=dat[,1],rep2=dat[,2],
                   rank1=rank(-dat[,1]),rank2=rank(-dat[,2]),
                   idr=res$idr)
  
  ggplot(df, aes(rep1,rep2,col=idr)) + geom_point()
  
  ggplot(df, aes(rank1,rank2,col=idr)) + geom_point()
  
}