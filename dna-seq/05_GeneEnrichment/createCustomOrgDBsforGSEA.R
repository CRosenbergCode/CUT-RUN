library(readr)
library("tools")
library(clusterProfiler)
library(tidyverse)
library(fgsea)
library(AnnotationForge)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


#The general workflow is as follows:
#Obtain GO annotations from Blast2GO, EGGNOG, or a similar platform
#Format GO annotations such that it is a two column file with 
#Split this file in smaller subsets that can be processed in parallel
#Use the GoArray.sh script to call () on these smaller files
#Concatenate the files with cat *. > GoAnnos.tsv
#Run the makeCustomOrgDB with your GoAnnos.tsv file and a csv annotation file containing names, symbols, and chromosomes

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
                 species="test",
                 goTable="go")
}

#Ex:
makeCustomOrgDB("Aedes_Ae_gene_annotations_PE20250124_GO2_14_2025.csv","EggnogGo_2025_02_21.tsv")


makeCustomOrgDB("Aedes_Ae_gene_annotations_PE20250124_GO2_14_2025.csv","AllGoSorted_317.tsv")




##### THIS FUNCTION IS CURRENTLY DEPRECATED AND THE PURPOSE IS ACCOMPLISHED IN BASH INSTEAD

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






#customFile is a csv file that follows the format of Dr. Corey Rosenberg's custom annotations
#Most importantly, there is a column named "Gene.ID" which contains the gene symbols to be matched against (ex: AAEL)
#outFile is a string containing the name of the file to which the output will be written
#goFile is a 
addGoAnnotationsToCustom = function(customFile, goFile, outFile){
  
  originalFrame=read.csv(customFile)
  goCols=read.table(goFile,sep="\t",header = FALSE, na.strings = "", fill = TRUE)
  
  
  originalFrame[ ,'Go.Terms'] = NA
  
  print(nrow(goCols))
  for(i in seq(nrow(goCols))){
    if(nrow(originalFrame[originalFrame$Gene.ID==goCols[i,1],])>0){
      originalFrame[originalFrame$Gene.ID==goCols[i,1],]$Go.Terms=goCols[i,2]
    }
    else{
      print(goCols[i,1])
    }
    #Add debugging in?
  }
  write.csv(originalFrame,outFile,row.names = FALSE)
  
}

#customFile is a string containing the name of the csv file to read the annotations from
#outFile is a string containing the name of the file to which the output will be written
#columns is a vector of two columns. The first is the gene ID and the second is the custom annotation(s) assigned to that gene. 
#By default, these are 2 and 4 to match the format of Dr. Corey Rosenberg's custom annotations
getCustomAnnotations = function(customFile,outFile,columns=c(2,4)){
  originalFrame=read.csv(customFile)
  newFrame=originalFrame[columns]
  newFrame=newFrame[!duplicated(newFrame), ]
  write_tsv(newFrame,outFile)
}