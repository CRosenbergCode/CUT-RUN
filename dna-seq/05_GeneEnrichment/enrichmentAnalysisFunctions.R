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

#The following library is a custom package which must be manually installed from source. 
#This can be performed with the following command, where "org.Aaegypti.eg.db" is the path to the org.Aaegypti.eg.db folder
#install.packages("org.Aaegypti.eg.db", repos=NULL,type="source")
#This folder is found in the Rosenberg NAS Informatics Resources -> Aedes_aegypti
library(org.Aaegypti.eg.db)
#If using a different orgDB than Aedes aegypti, it must be installed and loaded as seen above


#rnaFile is a csv which contains the results of differential gene expression. Generally the output of DEseq analysis
#rankMetric is the the choice of ranking for . Currently only supports a single, but more will be added in future
#goCat represents which of the three subontologies (CC,MF,BP) to include in the analysis. All three are included by default
rnaGSEA = function(rnaFile,orgData,rankMetric = "pval",goCat="All",minSize=15,maxSize=500,rankedList=c()){
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
  #return(newRank_pvalueAndFC)
  #myDEresults=read.csv(rnaFile)
  #myDEresults=myDEresults[!is.na(myDEresults$padj),]
  #newRank_pvalueAndFC = -log10(myDEresults$padj) *sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
  #names(newRank_pvalueAndFC) = myDEresults$geneID
  #newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
  set.seed(2025)
  egoCC <- gseGO(geneList     = newRank_pvalueAndFC,
                 OrgDb        = orgData,
                 keyType="GID",
                 ont          = goCat,
                 minGSSize    = minSize,
                 maxGSSize    = maxSize,
                 pvalueCutoff = 0.5,
                 verbose      = TRUE)
  
  head(egoCC)
  
  goplot(egoCC)
  dotplot(egoCC)+theme(axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"))#+theme(text = element_text(size = ))
  return(egoCC)
}

#Example Usage
exampleGo=rnaGSEA("RNAseqPilot.csv",goCat="BP",orgData = org.Aaegypti.eg.db)

head(exampleGo)

#Plotting 
goplot(exampleGo,showCategory = 5)
dotplot(exampleGo)
#Showing first few GO categories that are most significantly differentially expressed
head(exampleGo)


sfGo=rnaGSEA("DAY1_BFvSFdresultsp10_ms-PEannotated.csv",goCat="MF",orgData = org.Aaegypti.eg.db)
head(sfGo)
goplot(sfGo)
dotplot(sfGo)



#Using log2 fold change instead of p-value
fc=rnaGSEA("DAY1_BFvSFdresultsp10_ms-PEannotated.csv",goCat="MF",rankMetric ="log2fc",orgData = org.Aaegypti.eg.db)
goplot(fc)
dotplot(fc)
head(fc)

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

myGMT=gmtGSEA("RNAseqPilot.csv","RosenbergCustomAnnotations_2025_2_28.gmt")
head(myGMT)

#Return multiple objects used for plotting

multiGMT=gmtGSEA("RNAseqPilot.csv","RosenbergCustomAnnotations_2025_2_28.gmt",returnMultiple=TRUE)
#head(myGMT)

plotEnrichment(multiGMT$paths[["IMM"]],
               multiGMT$genes) + labs(title="Immune Genes")


#Enrichment Score Plot Subsetted by significance
multiGMT=gmtGSEA("RNAseqPilot.csv","RosenbergCustomAnnotations_2025_2_28.gmt",returnMultiple=TRUE)
selected_rows=multiGMT$table[multiGMT$table$padj<0.05,]
selected_rows = selected_rows %>%
  arrange(padj)
selected_GO=multiGMT$paths[selected_rows$pathway]
selected_rows$pathway

plotGseaTable(selected_GO, multiGMT$genes, selected_rows, 
              gseaParam=0.5,colwidths =  c(2, 5, 0.8, 0, 1.2))+theme(axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"))

#TODO\
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


testCompareRNA=rnaCompareGSEA(c("DAY1_BFvSFdresultsp10_ms-PEannotated.csv","RNAseqPilot.csv"),orgData=org.Aaegypti.eg.db)

dotplot(testCompareRNA,showCategory = 5)


#"DAY1_BFvSFdresultsp10_ms-PEannotated.csv"
#RNAseqPilot.csv

chipORA = function(peakFile,orgData,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",goCat="All",qval=0.5){
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
                  qvalueCutoff = qval, 
                  readable = TRUE)
  return(ego)
}

#Example Usage


testORA=chipORA("C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\BF_Ac_Rep1_Control_peaks.narrowPeak",org.Aaegypti.eg.db,goCat="BP")

testMF=chipORA("C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\BF_Ac_Rep1_Control_peaks.narrowPeak",org.Aaegypti.eg.db,goCat="MF")


head(testORA)

dotplot(testORA,showCategory = 5)
cnetplot(testORA,showCategory = 2)

head(testMF)
dotplot(testMF,showCategory = 5)
cnetplot(testMF,showCategory = 5)



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

compareChipGO=function(fileList,orgData,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",goCat="All",qval=0.5,plotting=FALSE,sampNames=c()){
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
  pre_genes = lapply(exampleDistances, function(i) as.data.frame(i)$geneId)

  names(pre_genes)=sampNames
  #genes = list(X1=unique(pre_genes[[1]]),X2=unique(pre_genes[[2]]))#list(X1=exampleDistances[1],X2=exampleDistances[2])
  compGO <- compareCluster(geneCluster = pre_genes,#genes, 
                           #universe=backgroundids,
                           fun = "enrichGO",
                           keyType = "GID", 
                           OrgDb = orgData, 
                           ont = goCat, 
                           qvalueCutoff = qval, 
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


testcompare2=compareChipGO(samplefiles2,org.Aaegypti.eg.db,goCat="MF")

dotplot(testcompare2,showCategory = 5)


testcompare5=compareChipGO(samplefiles5,org.Aaegypti.eg.db,goCat="MF",sampNames = c("BF_Ac_1","BF_Ac_2","RVFV_Ac_1_1","RVFV_1_1","RVFV_Ac_2"))

dotplot(testcompare5,showCategory = 5)

head(testcompare5)


metaFile=read.csv("CUT_RUN_Meta_File_MACS2_0.05_keepdup.csv")
d7_me_Pilot=metaFile[metaFile$RiftExperiment==TRUE,]
d7_me_Pilot=d7_me_Pilot[d7_me_Pilot$Factor=="H3K9Me3",]

samplefiles=d7_me_Pilot$Peaks


testcompare8=compareChipGO(samplefiles,org.Aaegypti.eg.db,goCat="MF")#,sampNames = c("BF_Ac_1","BF_Ac_2","RVFV_Ac_1_1","RVFV_1_1","RVFV_Ac_2"))

dotplot(testcompare8,showCategory = 5)


#Plotting to look at further

#p1 <- heatplot(edox, showCategory=5)
#p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
#cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
