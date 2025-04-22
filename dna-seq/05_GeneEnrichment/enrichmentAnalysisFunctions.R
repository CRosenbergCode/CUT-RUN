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

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.20")


#Necessary for Windows machines
register(SerialParam())


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#The following library is a custom package which must be manually installed from source. 
#This can be performed with the following command, where "org.Aaegypti.eg.db" is the path to the org.Aaegypti.eg.db folder
#install.packages("org.Aaegypti.eg.db", repos=NULL,type="source")
#This folder is found in the Rosenberg NAS Informatics Resources -> Aedes_aegypti
library(org.Aaegypti.eg.db)
#If using a different orgDB than Aedes aegypti, it must be installed and loaded as seen above


##FUNCTIONS IN THIS FILE

#rnaGSEA: perform GSEA analysis on single file, which should be differential expression output (generally DESeq2)
#Should be used for GO annotations
#Requires an OrgDB library

#gmtGSEA: As above, but using a GMT file instead of a GMT. Should be used for custom annotations
#Requires a GMT file

#rnaCompareGSEA: As with rnaGSEA, but with multiple files instead of one. Should be used to compare results side-by-side
#PLEASE NOTE THIS HELPS VISUALIZE MULTIPLE SAMPLES, BUT RESULTS ARE SIMPLY GSEA FOR INDIVIDUAL SAMPLES. THERE IS NO DIRECT STATISTICAL COMPARISON BETWEEN SAMPLES
#Requires an OrgDB library

#chipORA: Perform overrepresentation analysis using the presence of peaks, based on a peak file
#Requires a gff file that includes all genomic features, not just genes and an OrgDB library
#Will fail without ability to contact NCBI website through internet

#compareChipGO: As above, but with multiple peak files instead of one. Should be used to compare results side-by-side
#PLEASE NOTE THIS HELPS VISUALIZE MULTIPLE SAMPLES, BUT RESULTS ARE SIMPLY GSEA FOR INDIVIDUAL SAMPLES. THERE IS NO DIRECT STATISTICAL COMPARISON BETWEEN SAMPLES
#Requires a gff file that includes all genomic features, not just genes and an OrgDB library
#Will fail without ability to contact NCBI website through internet

#presenceAbsenceORA: Perform overrepresentation analysis using the set of genes which have peaks in one file but not another
#Unlike the previous function, this can be used to directly compare two samples
#Requires a gff file that includes all genomic features, not just genes and an OrgDB library
#Will fail without ability to contact NCBI website through internet

#chipGSEA: perform GSEA analysis on single file, which should be differential binding output (generally Diffbind)
#Should be used for GO annotations
#Requires a gff file that includes all genomic features, not just genes and an OrgDB library
#Will fail without ability to contact NCBI website through internet


#Generally, plotting functions should be useable on the results of any of these functions








#rnaFile is a csv which contains the results of differential gene expression. Generally the output of DEseq analysis
#rankMetric is the the choice of ranking for . Currently only supports a single, but more will be added in future
#goCat represents which of the three subontologies (CC,MF,BP) to include in the analysis. All three are included by default
rnaGSEA = function(rnaFile,orgData,rankMetric = "pval",goCat="All",minSize=15,maxSize=500,rankedList=c(),keycol="GID"){
  myDEresults=read.csv(rnaFile)
  if(rankMetric == "pval"){
    myDEresults=read.csv(rnaFile)
    myDEresults=myDEresults[!is.na(myDEresults$pvalue),]
    newRank_pvalueAndFC = -log10(myDEresults$pvalue) *sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
    names(newRank_pvalueAndFC) = myDEresults$geneID
    newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
  }
  if(rankMetric == "log2fc"){
    myDEresults=read.csv(rnaFile)
    myDEresults=myDEresults[!is.na(myDEresults$pvalue),]
    newRank_pvalueAndFC = myDEresults$log2FoldChange #*sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
    names(newRank_pvalueAndFC) = myDEresults$geneID
    newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
  }
  if(rankMetric == "both"){
    myDEresults=read.csv(rnaFile)
    myDEresults=myDEresults[!is.na(myDEresults$pvalue),]
    newRank_pvalueAndFC = myDEresults$log2FoldChange * -log10(myDEresults$pvalue) #*sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
    names(newRank_pvalueAndFC) = myDEresults$geneID
    newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
  }

  set.seed(2025)
  egoCC <- gseGO(geneList     = newRank_pvalueAndFC,
                 OrgDb        = orgData,
                 keyType=keycol,
                 ont          = goCat,
                 minGSSize    = minSize,
                 maxGSSize    = maxSize,
                 pvalueCutoff = 0.5,
                 verbose      = TRUE)
  
  
  #head(egoCC)
  
  #goplot(egoCC)
  #dotplot(egoCC)+theme(axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"))#+theme(text = element_text(size = ))
  return(egoCC)
}

library(org.Ctarsalis.eg.db)
exampleGo=rnaGSEA("CulexAdultPlaqueTable_4_2_2025.csv",goCat="MF",orgData = org.Ctarsalis.eg.db,rankMetric ="both")

examplePCR=rnaGSEA("CulexAdultPCRTable_4_2_2025.csv",goCat="MF",orgData = org.Ctarsalis.eg.db,rankMetric ="both")

dotplot(examplePCR)
testTemp=read.csv("CulexAdultPlaqueTable_4_2_2025.csv")
write.csv(testTemp[testTemp$pvalue<0.9,],"CulexAdultPCRTable_4_2_2025.5.csv",row.names = FALSE)

examplePCR=rnaGSEA("CulexAdultPCRTable_4_2_2025.5.csv",goCat="MF",orgData = org.Ctarsalis.eg.db,rankMetric ="both")

#Plotting 
goplot(exampleGo,showCategory = 5)
dotplot(exampleGo)

#CulexAdultPlaqueTable_4_2_2025.csv

#Example Usage
exampleGo=rnaGSEA("RNAseqPilot.csv",goCat="BP",orgData = org.Aaegypti.eg.db)


options(enrichplot.colours = c("yellow","purple"))
options(enrichplot.colours = c("red","blue"))
dotplot(exampleGo)

head(exampleGo)

#Plotting 
goplot(exampleGo,showCategory = 5)
dotplot(exampleGo)
#Showing first few GO categories that are most significantly differentially expressed
head(exampleGo)



#Using log2 fold change instead of p-value
fc=rnaGSEA("DAY1_BFvSFdresultsp10_ms-PEannotated.csv",goCat="MF",rankMetric ="log2fc",orgData = org.Aaegypti.eg.db)
goplot(fc)
dotplot(fc)
example=head(fc)

#Using the product of log2 fold change and -log10(pvalue)
bothVals=rnaGSEA("DAY1_BFvSFdresultsp10_ms-PEannotated.csv",goCat="MF",rankMetric ="both",orgData = org.Aaegypti.eg.db)
goplot(bothVals)
dotplot(bothVals)
head(bothVals)

#chipGSEA = function()

# preparing geneSet collections...
# GSEA analysis...
# no term enriched under specific pvalueCutoff...
# Error in .checkKeysAreWellFormed(keys) : 
#   keys must be supplied in a character vector with no NAs
# In addition: Warning message:
#   In preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
#                                There are ties in the preranked stats (15.69% of the list).
#                              The order of those tied genes will be arbitrary, which may produce unexpected results.
#                              Called from: .checkKeysAreWellFormed(keys)


#Add the no gene ID error



#TODO
#Write description

gmtGSEA = function(rnaFile,gmtFile,returnMultiple=FALSE){
  DE_results=read.csv(rnaFile)
  DE_results=DE_results[!is.na(DE_results$pvalue),]
  gene_list = -log10(DE_results$pvalue) *sign(DE_results$log2FoldChange)
  names(gene_list) = DE_results$geneID
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  head(gene_list)
  gene_list
  sum(is.na(gene_list))
  
  
  myGO = fgsea::gmtPathways(gmtFile)
  set.seed(2025)
  fgResPVal <- fgsea::fgsea(pathways =myGO, 
                            stats = gene_list,
                            nPermSimple = 10000,
                            eps = 0) 
  writeP <- apply(fgResPVal,2,as.character)
  if(returnMultiple){
    return(list(table = fgResPVal,paths=myGO,genes=gene_list))
  }
  return(fgResPVal)

}

#Example Usage

#Return just table of results

myGMT=gmtGSEA("RVFVvBF_DAY1_DESeqresults_min_10_all-annot.csv","AAE_2025-4-16.gmt",returnMultiple = FALSE)

head(myGMT)
goplot(myGMT)

#Return multiple objects used for plotting

multiGMT=gmtGSEA("RNAseqPilot.csv","RosenbergCustomAnnotations_2025_3_24.gmt",returnMultiple=TRUE)
#head(myGMT)

plotEnrichment(multiGMT$paths[["IMM"]],
               multiGMT$genes) + labs(title="Immune Genes")


#Enrichment Score Plot Subsetted by significance
multiGMT=gmtGSEA("RNAseqPilot.csv","RosenbergCustomAnnotations_2025_3_24.gmt",returnMultiple=TRUE)


head(multiGMT$table)

selected_rows=multiGMT$table[multiGMT$table$padj<0.05,]
selected_rows = selected_rows %>%
  arrange(padj)
selected_GO=multiGMT$paths[selected_rows$pathway]

plotGseaTable(selected_GO, multiGMT$genes, selected_rows, 
              gseaParam=0.5,colwidths =  c(2, 5, 0.8, 0, 1.2))+theme(axis.text.y = element_text(color = "grey20", size = 50, angle = 0, hjust = 1, vjust = 0, face = "plain"))
#Size was 15
#TODO
#Write description


rnaCompareGSEA = function(fileList,orgData,rankMetric = "pval",goCat="All",pval=0.5,plotting=FALSE,sampNames=c(),minSize=15,maxSize=500){
  getRankedGenes = function(rnaFile){
    myDEresults=read.csv(rnaFile)
    if(rankMetric == "pval"){
      myDEresults=read.csv(rnaFile)
      myDEresults=myDEresults[!is.na(myDEresults$pvalue),]
      newRank_pvalueAndFC = -log10(myDEresults$pvalue) *sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
      names(newRank_pvalueAndFC) = myDEresults$geneID
      newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
    }
    if(rankMetric == "log2fc"){
      myDEresults=read.csv(rnaFile)
      myDEresults=myDEresults[!is.na(myDEresults$pvalue),]
      newRank_pvalueAndFC = myDEresults$log2FoldChange #*sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
      names(newRank_pvalueAndFC) = myDEresults$geneID
      newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
    }
    if(rankMetric == "both"){
      myDEresults=read.csv(rnaFile)
      myDEresults=myDEresults[!is.na(myDEresults$pvalue),]
      newRank_pvalueAndFC = myDEresults$log2FoldChange * -log10(myDEresults$pvalue) #*sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
      names(newRank_pvalueAndFC) = myDEresults$geneID
      newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
    }
    

    return(newRank_pvalueAndFC)
  }
  
  
  
  
  if(length(sampNames)==0){
    sampNames=paste0("Sample_", 1:length(fileList))
  }
  
  pre_genes = lapply(fileList,getRankedGenes)#function(i) as.data.frame(i)$geneId)
  names(pre_genes)=sampNames
  compGO <- compareCluster(geneCluster = pre_genes,
                           fun = "gseGO",
                           keyType = "GID", 
                           OrgDb = orgData, 
                           ont = goCat, 
                           minGSSize    = minSize,
                           maxGSSize    = maxSize,
                           pvalueCutoff = pval,
                           pAdjustMethod = "BH")#, readable=TRUE)
  if(plotting){
    dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis") 
  }
  return(compGO)
  return 
}


#Example Usage
testCompareRNA=rnaCompareGSEA(c("DAY1_BFvSFdresultsp10_ms-PEannotated.csv","RNAseqPilot.csv"),orgData=org.Aaegypti.eg.db)

dotplot(testCompareRNA,showCategory = 5)

testCompareRNA@compareClusterResult


#https://github.com/alserglab/fgsea/issues/29

#TODO
#Write description

chipORA = function(peakFile,orgData,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",goCat="All",qval=0.5,promoterOnly=TRUE){
  plotPeakDistances = function(peakFileList,sampNames=c(),verbose=TRUE,txFile=txFile){
    samplefiles <- as.list(peakFileList)
    
    testtx=makeTxDbFromGFF(txFile,format="gff3")
    
    
    peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=testtx, 
                           tssRegion=c(-2000, 2000), verbose=FALSE)
    
    
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
  
  exampleDistances=plotPeakDistances(c(peakFile),txFile=txFile)
  
  ex_annot = exampleDistances[[1]]@anno
  
  if(promoterOnly){
    ex_annot=ex_annot[abs(ex_annot$distanceToTSS) <= 2000,]
  }
  
  gIDs <- ex_annot$geneId
  
  annotations_edb <- AnnotationDbi::select(orgData,
                                           keys = gIDs,
                                           columns = c("GENENAME"))#,
  
  annotations_edb$GID <- as.character(annotations_edb$GID)
  
  
  
  ego <- enrichGO(gene = gIDs, 
                  keyType = "GID", 
                  OrgDb = orgData, 
                  ont = goCat, 
                  pAdjustMethod = "BH", 
                  pvalueCutoff = qval,
                  readable = TRUE)
  return(ego)
}

#Example Usage


testORA=chipORA("C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\BF_Ac_Rep1_Control_peaks.narrowPeak",org.Aaegypti.eg.db,goCat="BP")
dotplot(testORA,showCategory = 5)


testORA=chipORA("C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\BF_Ac_Rep1_Control_peaks.narrowPeak",org.Aaegypti.eg.db,goCat="BP",promoterOnly=FALSE)
dotplot(testORA,showCategory = 5)


testMF=chipORA("C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\BF_Ac_Rep1_Control_peaks.narrowPeak",org.Aaegypti.eg.db,goCat="MF",promoterOnly = FALSE)


head(testORA)

dotplot(testORA,showCategory = 5)
cnetplot(testORA,showCategory = 2)

head(testMF)
dotplot(testMF,showCategory = 5)
cnetplot(testMF,showCategory = 5)

testMF=chipORA("C:\\Users\\hunte\\Desktop\\AltChip\\PilotAcBF50K_200K_Distal.peakFile",org.Aaegypti.eg.db,goCat="MF",promoterOnly = FALSE)


#Has very small text by default, will need to manually adjust size per plot
#Size of test is inversely proportional to the value of cex. Return to default of 1 after calling
par(cex = 0.25)
plotGOgraph(testORA)
par(cex = 1)
goplot(testORA)



#fileList is a vector with each element being the path to a 
#rankMetric is the the choice of ranking for . Currently only supports a single, but more will be added in future
#goCat represents which of the three subontologies (CC,MF,BP) to include in the analysis. All three are included by default
#orgData is an organism db object, which should be referred to using the same name as its library

compareChipGO=function(fileList,orgData,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",goCat="All",qval=0.5,plotting=FALSE,sampNames=c(),promoterOnly=TRUE){
  plotPeakDistances = function(peakFileList,sampNames=c(),verbose=TRUE,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff"){
    samplefiles <- as.list(peakFileList)
    
    testtx=makeTxDbFromGFF(txFile,format="gff3")
    
    
    peakAnnoList <- lapply(peakFileList, annotatePeak, TxDb=testtx, 
                           tssRegion=c(-2000, 2000), verbose=FALSE)
    
    
    if(verbose){
      print(peakAnnoList)
    }
    if(length(sampNames)>0){
      names(peakAnnoList)=sampNames
    }
    par("mar")
    par(mar=c(1,1,1,1))
    plotAnnoBar(peakAnnoList)
    plotDistToTSS(peakAnnoList, title="Distribution of Peaks\n relative to TSS")
    return(peakAnnoList)
  }
  if(length(sampNames)==0){
    sampNames=paste0("Sample_", 1:length(fileList))
  }
  exampleDistances=plotPeakDistances(fileList,txFile=txFile,sampNames = sampNames)
  
  if(promoterOnly){
    pre_genes = lapply(exampleDistances, function(i) as.data.frame(i))
    for(i in seq(length(pre_genes))){
      head(pre_genes[[i]])
      pre_genes[[i]]=pre_genes[[i]][abs(pre_genes[[i]]$distanceToTSS) <= 2000,]
    }
    pre_genes = lapply(pre_genes, function(i) as.data.frame(i)$geneId)
  }
  else{
    pre_genes = lapply(exampleDistances, function(i) as.data.frame(i)$geneId)
  }

  # if(promoterOnly){
  #   for(i in seq(length(pre_genes))){
  #     head(pre_genes[[i]])
  #     pre_genes[[i]]=pre_genes[[i]][abs(pre_genes[[i]]$distanceToTSS) <= 2000,]
  #   }
  # }
  
  names(pre_genes)=sampNames
  #genes = list(X1=unique(pre_genes[[1]]),X2=unique(pre_genes[[2]]))#list(X1=exampleDistances[1],X2=exampleDistances[2])
  compGO <- compareCluster(geneCluster = pre_genes,#genes, 
                           #universe=backgroundids,
                           fun = "enrichGO",
                           keyType = "GID", 
                           OrgDb = orgData, 
                           ont = goCat, 
                           pvalueCutoff = qval,
                           pAdjustMethod = "BH")#, readable=TRUE)
  if(plotting){
    dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis") 
  }
  return(compGO)
}

#Example Usage

samplefiles2=c("C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\BF_Ac_Rep1_Control_peaks.narrowPeak","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\RVFV_Ac_Rep1_1_Control_peaks.narrowPeak")#RVFV_Ac_Rep1_2_Control_peaks.narrowPeak")#,"C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewSEACRPeaks\\RVFV_D7_Ac_Rep1_1.relaxed.bed","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewSEACRPeaks\\RVFV_D7_Ac_Rep1_2.relaxed.bed","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewSEACRPeaks\\RVFV_D7_Ac_Rep2.relaxed.bed")

samplefiles5=c("C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\BF_Ac_Rep1_Control_peaks.narrowPeak","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\BF_Ac_Rep2_Control_peaks.narrowPeak",
              "C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\RVFV_Ac_Rep1_1_Control_peaks.narrowPeak","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\RVFV_Ac_Rep1_2_Control_peaks.narrowPeak",
              "C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\RVFV_Ac_Rep2_Control_peaks.narrowPeak")

#Compare two files
testcompare2=compareChipGO(samplefiles2,org.Aaegypti.eg.db,goCat="MF")

dotplot(testcompare2,showCategory = 5)

#Compare five files
testcompare5=compareChipGO(samplefiles5,org.Aaegypti.eg.db,goCat="MF",sampNames = c("BF_Ac_1","BF_Ac_2","RVFV_Ac_1_1","RVFV_1_1","RVFV_Ac_2"))

dotplot(testcompare5,showCategory = 5)

head(testcompare5)

#Use the metadata file to grab peaks data for 8 files
metaFile=read.csv("CUT_RUN_Meta_File_MACS2_0.05_keepdup.csv")
d7_me_Pilot=metaFile[metaFile$RiftExperiment==TRUE,]
d7_me_Pilot=d7_me_Pilot[d7_me_Pilot$Factor=="H3K9Me3",]


testcompare8=compareChipGO(samplefiles,org.Aaegypti.eg.db,goCat="MF")#,sampNames = c("BF_Ac_1","BF_Ac_2","RVFV_Ac_1_1","RVFV_1_1","RVFV_Ac_2"))

dotplot(testcompare8,showCategory = 5)

testcompareDistalPilot=compareChipGO(c("PilotAcRVFV50K_200K_Distal.peakFile","PilotAcBF50K_200K_Distal.peakFile"),org.Aaegypti.eg.db,goCat="MF",promoterOnly = FALSE,sampNames = c("Distal RVFV_Ac_D7","Distal BF_Ac_D7"))

dotplot(testcompareDistalPilot,showCategory = 5)


#TODO
#Write description

#Implemented, but not inspiring results

presenceAbsenceORA = function(peakFile1,peakFile2,orgData,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",goCat="All",qval=0.5,plotting=FALSE,sampNames=c(),promoterOnly=TRUE){
  
  plotPeakDistances = function(peakFileList,sampNames=c(),verbose=TRUE,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff"){
    samplefiles <- as.list(peakFileList)
    
    testtx=makeTxDbFromGFF(txFile,format="gff3")
    
    
    peakAnnoList <- lapply(peakFileList, annotatePeak, TxDb=testtx, 
                           tssRegion=c(-2000, 2000), verbose=FALSE)
    
    
    if(verbose){
      print(peakAnnoList)
    }
    if(length(sampNames)>0){
      names(peakAnnoList)=sampNames
    }
    par("mar")
    par(mar=c(1,1,1,1))
    plotAnnoBar(peakAnnoList)
    plotDistToTSS(peakAnnoList, title="Distribution of Peaks\n relative to TSS")
    return(peakAnnoList)
  }
  if(length(sampNames)==0){
    sampNames=paste0("Sample_", 1:2)
  }
  samp1=as.data.frame(plotPeakDistances(peakFile1,txFile=txFile)[[1]])
  samp2=as.data.frame(plotPeakDistances(peakFile2,txFile=txFile)[[1]])
  if(promoterOnly){
    samp1_peaks=samp1[abs(samp1$distanceToTSS) <= 2000,]$geneId
    samp2_peaks=samp2[abs(samp2$distanceToTSS) <= 2000,]$geneId
  }
  else{
    samp1_peaks=samp1$geneId
    samp2_peaks=samp2$geneId
  }
  #samp1_peaks=samp1[abs(samp1$distanceToTSS) <= 2000,]$geneId
  
  #samp2=as.data.frame(testcompare2$Sample_2)
  #samp2_peaks=samp2[abs(samp2$distanceToTSS) <= 2000,]$geneId
  diffPeaks=samp1_peaks[!samp1_peaks %in% samp2_peaks]
  #return(diffPeaks)
  
  ego <- enrichGO(gene = diffPeaks, 
                  keyType = "GID", 
                  OrgDb = orgData, 
                  ont = goCat, 
                  pAdjustMethod = "BH", 
                  pvalueCutoff = qval,
                  #qvalueCutoff = qval, 
                  readable = TRUE)
  return(ego)
}


#Example usage
metaFile=read.csv("CUT_RUN_Meta_File_MACS2_0.05_keepdup.csv")
d7_me_Pilot=metaFile[metaFile$RiftExperiment==TRUE,]
d7_me_Pilot=d7_me_Pilot[d7_me_Pilot$Factor=="H3K9Me3",]
d7_me_Pilot_Peaks=d7_me_Pilot$Peaks

metaFile=read.csv("CUT_RUN_Meta_File_MACS2_0.05_keepdup.csv")
d7_ac_Pilot=metaFile[metaFile$RiftExperiment==TRUE,]
d7_ac_Pilot=d7_ac_Pilot[d7_ac_Pilot$Factor=="H3K27Ac",]
d7_ac_Pilot_Peaks=d7_ac_Pilot$Peaks



testPresenceAbsence=presenceAbsenceORA(d7_me_Pilot_Peaks[3],d7_me_Pilot_Peaks[4],orgData=org.Aaegypti.eg.db,goCat="MF",promoterOnly=FALSE)

testPresenceAbsence=presenceAbsenceORA(d7_me_Pilot_Peaks[3],d7_me_Pilot_Peaks[4],orgData=org.Aaegypti.eg.db,goCat="MF",promoterOnly=FALSE)


dotplot(testPresenceAbsence)

head(testPresenceAbsence)


#TODO
#Write documentation
ChIPGSEA=function(peakFile,orgData,rankMetric = "pval",txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",goCat="All",qval=0.5,promoterOnly=TRUE,enhancerOnly=FALSE,minSize=15,maxSize=500,plotting=TRUE,dists=c(-2000,2000)){
  
  tempFile="tempDiffbind.tsv"
  
  testtx=makeTxDbFromGFF(txFile)
  
  myDiffbind=read.csv(peakFile)
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
  print("The number of rows after filtering is")
  print(nrow(ex_annot))
  #Initially sorted by -log10(pvalue) to keep most significant gene for category, but could change for other methods 
  ex_annot$pValue
  ex_annot=ex_annot %>%
    arrange(desc(pValue))
  
  unique_annot=ex_annot[!duplicated(ex_annot$geneId),]
  
  print("After adjusting for uniqueness")
  print(nrow(unique_annot))
  
  if(rankMetric=="pval"){
    unique_annot=arrange(unique_annot,desc(unique_annot$pValue*sign(unique_annot$signalValue)))
    
    myranked=unique_annot$pValue*sign(unique_annot$signalValue)
  }
  
  if(rankMetric=="foldchange"){
    unique_annot=arrange(unique_annot,desc(unique_annot$signalValue))
    
    myranked=unique_annot$signalValue
  }
  
  if(rankMetric=="both"){
    unique_annot=arrange(unique_annot,desc(unique_annot$pValue*unique_annot$signalValue))
    myranked=unique_annot$pValue*unique_annot$signalValue
  }
  print("Ready to GO!")
  
  names(myranked)=unique_annot$geneId
  
  #Significantly better results with duplicate peaks
  set.seed(2025)
  ego <- gseGO(geneList     = myranked,
               OrgDb        = orgData,
               keyType="GID",
               ont          = goCat,
               minGSSize    = minSize,
               maxGSSize    = maxSize,
               pvalueCutoff = 0.9,
               verbose      = TRUE)
  if(plotting){
    dotplot(ego) 
  }
  
  if(file.exists(tempFile)) {
    file.remove(tempFile)
  }
  return(ego)
}

#Example Usage
#The input file is the output of diffbind, written to a csv

#Using only promoter regions
mychip=ChIPGSEA("PilotNewMerged_Results.csv",orgData = org.Aaegypti.eg.db,goCat = "BP")


#Using only putative enhancer regions
bpchip=ChIPGSEA("PilotNewMerged_Results.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)

dotplot(bpchip)

bpchip=ChIPGSEA("BF_Ac_D1vD3DiffbindResults.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = TRUE,enhancerOnly = FALSE)

dotplot(bpchip)

compareChIPGSEA=function(fileList,orgData,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",rankMetric = "pval",goCat="All",qval=0.5,plotting=FALSE,sampNames=c(),promoterOnly=TRUE,enhancerOnly=FALSE,dists=c(-2000,2000),keycol="GID"){
  #myRankedGenes=vector("list", length = length(fileList))
  tempFile="tempDiffbind.tsv"
  testtx=makeTxDbFromGFF(txFile)
  count=0
  chipgseahelper=function(diffbindFile){
    
    myDiffbind=read.csv(diffbindFile)
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
    print("The number of rows after filtering is")
    print(nrow(ex_annot))
    #Initially sorted by -log10(pvalue) to keep most significant gene for category, but could change for other methods 
    ex_annot$pValue
    ex_annot=ex_annot %>%
      arrange(desc(pValue))
    
    unique_annot=ex_annot[!duplicated(ex_annot$geneId),]
    
    print("After adjusting for uniqueness")
    print(nrow(unique_annot))
    
    if(rankMetric=="pval"){
      unique_annot=arrange(unique_annot,desc(unique_annot$pValue*sign(unique_annot$signalValue)))
      
      myranked=unique_annot$pValue*sign(unique_annot$signalValue)
    }
    
    if(rankMetric=="foldchange"){
      unique_annot=arrange(unique_annot,desc(unique_annot$signalValue))
      
      myranked=unique_annot$signalValue
    }
    
    if(rankMetric=="both"){
      unique_annot=arrange(unique_annot,desc(unique_annot$pValue*unique_annot$signalValue))
      myranked=unique_annot$pValue*unique_annot$signalValue
    }
    names(myranked)=unique_annot$geneId
    return(myranked)
    #myRankedGenes[count]=myranked
    #count=count+1
  }
  
  #names(my)
  #Significantly better results with duplicate peaks
  # ego <- gseGO(geneList     = myRankedGenes,
  #              OrgDb        = orgData,
  #              keyType="GID",
  #              ont          = goCat,
  #              minGSSize    = minSize,
  #              maxGSSize    = maxSize,
  #              pvalueCutoff = 0.9,
  #              verbose      = TRUE)
  #myRankedGenes=list(myRankedGenes)
  #names(myRankedGenes)=sampNames
  #print(myRankedGenes)
  pre_genes = lapply(fileList,chipgseahelper)#function(i) as.data.frame(i)$geneId)
  
  if(length(sampNames)==0){
    sampNames=paste0("Sample_", 1:length(fileList))
  }
  
  names(pre_genes)=sampNames
  set.seed(2025)
  compGO <- compareCluster(geneCluster = pre_genes,#genes, 
                           #universe=backgroundids,
                           fun = "gseGO",
                           keyType = "GID", 
                           OrgDb = orgData, 
                           ont = goCat, 
                           pvalueCutoff = qval,
                           pAdjustMethod = "BH")#, readable=TRUE)
  if(plotting){
    dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis") 
  }
  if(file.exists(tempFile)) {
    file.remove(tempFile)
  }
  return(compGO)
}

#Example usage
compareResults=compareChIPGSEA(c("BF_Ac_D1vD3DiffbindResults.csv","BF_Ac_D3vD7DiffbindResults.csv"),orgData=org.Aaegypti.eg.db,goCat="BP",sampNames=c("BF D1vd3","BF D3vD7"))
dotplot(compareResults)


chipGMT = function(peakFile,gmtFile,rankMetric = "pval",txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",qval=0.5,promoterOnly=TRUE,dists=c(-2000,2000),returnMultiple=FALSE){
  tempFile="tempDiffbind.tsv"
  
  testtx=makeTxDbFromGFF(txFile)
  
  myDiffbind=read.csv(peakFile)
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
  print("The number of rows after filtering is")
  print(nrow(ex_annot))
  #Initially sorted by -log10(pvalue) to keep most significant gene for category, but could change for other methods 
  ex_annot$pValue
  ex_annot=ex_annot %>%
    arrange(desc(pValue))
  
  unique_annot=ex_annot[!duplicated(ex_annot$geneId),]
  
  print("After adjusting for uniqueness")
  print(nrow(unique_annot))

  print("Ready to GO!")
  
  if(rankMetric=="pval"){
    unique_annot=arrange(unique_annot,desc(unique_annot$pValue*sign(unique_annot$signalValue)))
    
    myranked=unique_annot$pValue*sign(unique_annot$signalValue)
  }
  
  if(rankMetric=="foldchange"){
    unique_annot=arrange(unique_annot,desc(unique_annot$signalValue))
    
    myranked=unique_annot$signalValue
  }
  
  if(rankMetric=="both"){
    unique_annot=arrange(unique_annot,desc(unique_annot$pValue*unique_annot$signalValue))
    myranked=unique_annot$pValue*unique_annot$signalValue
  }
  names(myranked)=unique_annot$geneId
  
  myGO = fgsea::gmtPathways(gmtFile)
  set.seed(2025)
  fgResPVal <- fgsea::fgsea(pathways =myGO, 
                            stats = myranked,
                            nPermSimple = 10000,
                            eps = 0) 
  writeP <- apply(fgResPVal,2,as.character)
  if(file.exists(tempFile)) {
    file.remove(tempFile)
  }
  if(returnMultiple){
    return(list(table = fgResPVal,paths=myGO,genes=myranked))
  }
  return(fgResPVal)
}

mychip=chipGMT("PilotNewMerged_Results.csv",gmtFile = "AAE_2025-4-19-9-50-9.gmt",returnMultiple = TRUE)

dotplot(mychip)

gmttable=mychip$table

selected_rows=mychip$table[mychip$table$padj<0.1,]
selected_rows = selected_rows %>%
  arrange(padj)
selected_GO=mychip$paths[selected_rows$pathway]
selected_rows$pathway

plotGseaTable(selected_GO, mychip$genes, selected_rows, 
              gseaParam=0.5,colwidths =  c(2, 5, 0.8, 0, 1.2))+theme(axis.text.y = element_text(color = "grey20", size = 50, angle = 0, hjust = 1, vjust = 0, face = "plain"))


#AAE_2025-4-16.gmt

rnaORA = function(rnaFile,orgData,goCat="All",minSize=15,maxSize=500,rankedList=c(),keycol="GID",asTable=FALSE,qval=0.5){
  if(!asTable){
    myDEresults=read.csv(rnaFile)
  }
  else{
    myDEresults=rnaFile
  }
  gIDs=myDEresults$geneID
  
  ego <- enrichGO(gene = gIDs, 
                keyType = "GID", 
                OrgDb = orgData, 
                ont = goCat, 
                pAdjustMethod = "BH", 
                pvalueCutoff = qval,
                readable = TRUE)
  return(ego)
}

pilot=read.csv("RNAseqPilot.csv")
pcr=read.csv("CulexAdultPCRTable_4_2_2025.5.csv")
subsetpcr=pcr[which(pcr$pvalue<0.005),]
subsetpilot=pilot[which(pilot$pvalue<0.01),]
egData=read.csv("Mock22_Mock30_DE_PVal.csv")
subseteg=egData[which(egData$padj<0.01)&which(egData$log2FoldChange<0),]
subseteg=egData[(egData$pvalue<0.005)&(egData$log2FoldChange>0)&(!is.na((egData$pvalue<0.005)&(egData$log2FoldChange>0))),]


pcrORA=rnaORA(subsetpcr,orgData = org.Ctarsalis.eg.db,asTable = TRUE,qval = 0.9)
head(pcrORA)

pcrORA=rnaORA(subsetpcr,orgData = org.Ctarsalis.eg.db,asTable = TRUE)

pilotORA=rnaORA(subsetpilot,orgData = org.Aaegypti.eg.db,asTable = TRUE)
head(pilotORA)


egORA=rnaORA(subseteg,orgData = org.Ctarsalis.eg.db,asTable = TRUE,qval = 0.9)
head(egORA)


