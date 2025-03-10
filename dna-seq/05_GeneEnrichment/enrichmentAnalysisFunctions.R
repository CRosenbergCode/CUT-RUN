library(readr)
library(clusterProfiler)
library(tidyverse)
library(fgsea)
library(stringr)
library(AnnotationForge)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

#The following library is a custom package which must be manually installed from source. 
#This can be performed with the following command, where "org.Aaegypti.eg.db" is the path to the org.Aaegypti.eg.db folder
#install.packages("org.Aaegypti.eg.db", repos=NULL,type="source")
#This folder is found in the Rosenberg NAS Informatics Resources -> Aedes_aegypti
library(org.Aaegypti.eg.db)


#rnaFile is a csv which contains the results of differential gene expression. Generally the output of DEseq analysis
#rankMetric is the the choice of ranking for . Currently only supports a single, but more will be added in future
#goCat represents which of the three subontologies (CC,MF,BP) to include in the analysis. All three are included by default
rnaGSEA = function(rnaFile,rankMetric = "padj",goCat="All",minSize=15,maxSize=500,rankedList=c()){
  myDEresults=read.csv(rnaFile)
  if(rankMetric == "padj"){
    myDEresults=read.csv(rnaFile)
    myDEresults=myDEresults[!is.na(myDEresults$padj),]
    newRank_pvalueAndFC = -log10(myDEresults$padj) *sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
    names(newRank_pvalueAndFC) = myDEresults$geneID
    newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
  }
  #myDEresults=read.csv(rnaFile)
  #myDEresults=myDEresults[!is.na(myDEresults$padj),]
  #newRank_pvalueAndFC = -log10(myDEresults$padj) *sign(myDEresults$log2FoldChange)#* abs(myDEresults$log2FoldChange)#
  #names(newRank_pvalueAndFC) = myDEresults$geneID
  #newRank_pvalueAndFC = newRank_pvalueAndFC[order(newRank_pvalueAndFC,decreasing = TRUE)]
  set.seed(2025)
  egoCC <- gseGO(geneList     = newRank_pvalueAndFC,
                 OrgDb        = org.Aaegypti.eg.db,
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
exampleGo=rnaGSEA("RNAseqPilot.csv",goCat="CC")

#Plotting 
goplot(exampleGo)
dotplot(exampleGo)
#Showing first few GO categories that are most significantly differentially expressed
head(exampleGo)


sfGo=rnaGSEA("DAY1_BFvSFdresultsp10_ms-PEannotated.csv",goCat="MF")
goplot(sfGo)
dotplot(sfGo)

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

gmtGSEA = function(rnaFile,gmtFile){
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
  fgResPVal <- fgsea::fgsea(pathways =myGO, 
                            stats = gene_list,
                            nPermSimple = 10000,
                            eps = 0) 
  writeP <- apply(fgResPVal,2,as.character)
  
  return(fgResPVal)
  #if(onlySig){
    
  #}

}