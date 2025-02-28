library(readr)
library(clusterProfiler)
library(tidyverse)
library(fgsea)
library(stringr)
library(AnnotationForge)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#The following library is a custom package which must be manually installed from source. 
#This can be performed with the following command, where "org.Aaegypti.eg.db" is the path to the org.Aaegypti.eg.db folder
#install.packages("org.Aaegypti.eg.db", repos=NULL,type="source")
library(org.Aaegypti.eg.db)


#rnaFile is a csv which contains the results of differential gene expression. Generally the output of DEseq analysis
#rankMetric is the the choice of ranking for . Currently only supports a single, but more will be added in future
#goCat represents which of the three subontologies (CC,MF,BP) to include in the analysis. All three are included by default
rnaGSEA = function(rnaFile,rankMetric = "default",goCat="All",minSize=15,maxSize=500){
  myDEresults=read.csv(rnaFile)
  if(rankMetric == "default"){
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

exampleGo=rnaGSEA("RNAseqPilot.csv",goCat="CC")
goplot(exampleGo)
dotplot(exampleGo)
sfGo=rnaGSEA("DAY1_BFvSFdresultsp10_ms-PEannotated.csv",goCat="CC")
goplot(sfGo)
dotplot(sfGo)

#chipGSEA = function()




#annotFile is a csv that contains gene names, locations, and more. This is based on the format of Dr. Rosenberg's custom annotations file
#goFile is a tsv with the first column of each row containing a single gene name (ex: AAEL00001) and the second containing a single go term (ex: GO:0000228)

makeCustomOrgDB = function(annotFile,goFile){
  annotFrame=read.csv(annotFile)
  
  annot.first <- annotFrame[match(unique(annotFrame$Gene.ID), annotFrame$Gene.ID),]
  mySym=annot.first[,c(2,8,7)]
  
  chrsCol=rep(NA,length(annot.first$Genomic.Location..Gene.))
  for(i in seq(length(annot.first$Genomic.Location..Gene.))){
    temp=unlist(strsplit(annot.first$Genomic.Location..Gene.[i],split=":"))[1]
    chrsCol[i]=temp
  }
  
  myChr=data.frame(Gene.ID=annot.first$Gene.ID,Chromosomes=chrsCol)
  myGo=read.table(goFile,header = TRUE)
  myEvidence=rep("IEA",length(myGo$Gene.ID))
  myGo["Evidence"]=myEvidence
  
  mySym=mySym[!duplicated(mySym),]
  myChr=myChr[!duplicated(myChr),]
  myGo=myGo[!duplicated(myGo),]
  names(mySym)=c("GID","SYMBOL","GENENAME")
  names(myChr)=c("GID","CHROMOSOME")
  names(myGo)=c("GID","GO","EVIDENCE")
  
  testGo=myGo[!is.na(myGo$GO),]
  
  testSym=mySym
  for(i in seq(length(testSym$SYMBOL))){
    if(testSym$SYMBOL[i]=="N/A"){
      testSym$SYMBOL[i]=testSym$GID[i]
    }
  }
  
  testIDs=testGo$GO
  
  count=0
  for(i in seq(length(testIDs))){
    if(!str_detect(testIDs[i], regex("GO:[0-9]+"))){
      testGo=testGo[-(i+count),]
      count=count+1
    }
  }
  
  ## Then call the function
  makeOrgPackage(gene_info=testSym, chromosome=myChr, go=testGo,
                 version="0.2",
                 maintainer="Hunter Ogg <hunter.a.ogg@gmail.com>",
                 author="Hunter Ogg <hunter.a.ogg@gmail.com>",
                 outputDir = "OrgDB",
                 tax_id="7159",
                 genus="Aedes",
                 species="aegypti",
                 goTable="go")
}

#Ex:
makeCustomOrgDB("Aedes_Ae_gene_annotations_PE20250124_GO2_14_2025.csv","EggnogGo_2025_02_21.tsv")

#The purpose of this function is to split tens of thousands of genes worth of annotations into smaller chunks that can be parallelized with GoArray.sh and GoArray.r

#file is a tsv with each row containing a single gene name and every go annotation corresponding to that gene. This is based on the format of eggnog output
#nGenes per file is the number of genes that will be written into each subfile. 
#base is the prefix of all of the files. It is "GoTermsSegment" by default
#

splitGoAnnos = function(file,nGenesPerFile=1000,base="GoTermsSegment"){
  annotFrame=read_tsv(file,col_names = c("Gene.ID","Go.Terms"))
  nSegments = floor(nrow(annotFrame)/nGenesPerFile)
  for(i in seq(nSegments)){
    first=(i-1)*1000+1
    last=i*1000
    subframe=annotFrame[first:last,]
    write.table(subframe,paste(base,i,".tsv",sep=""),sep="\t",row.names = FALSE,quote=FALSE)
  }
  testdf=annotFrame[(nSegments*nGenesPerFile+1):nrow(annotFrame),]
  print(nrow(testdf))
  write.table(testdf,paste(base,(nSegments+1),".tsv",sep=""),sep="\t",row.names = FALSE,quote=FALSE)
}

#Ex:
splitGoAnnos("MiniEggNogAedesAegypti68.tsv")