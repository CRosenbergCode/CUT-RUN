library(DESeq2) #Used as DESeq2 objects may be read in, although no library commands are invoked
library(readr) #For read_tsv()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#wholeGTF=read_tsv("AedesGenes.gtf",col_names = FALSE)

#Primarily a helper function which returns a list of the AAEL##### names from the gtf/gff file
#This is built for the Aedes aegypti vectorbase 10 column gff/gtf file specifically and may not work
#If applied to other organisms/formats as it relies upon a specific string format
getGeneName = function(gtfFrame,bedformat=FALSE){
  nsize=nrow(gtfFrame)
  retArr=rep("",nsize)
  for(i in seq(nsize)){
    retArr[i]=substring(toString(gtfFrame[i,10]),4,13)#,10,19)#,4,13)
  }
  return(retArr)
}

#The bulk of the code, this will be called by other functions after they have taken the information they need based
#on the type of input provided

#Can also be used directly with a vector of gene names
#Inputs
#Requires a gtf file, which is in the 10 column vectorbase format
#Requires a character vector, where each element is a 6 digital Aedes gene ID number (ex: AAEL009762)
#Returns a gtf with only the rows corresponding to the gene annotations provided, whether by file, deseq object, or vector
#If savename is provided, the resulting gtf will be saved with that name to the current working directory
subsetGeneNameHelper = function(gtfFrame,namesVec,savename=""){
  nsize=nrow(gtfFrame)
  subsetArr=rep(FALSE,nsize)
  geneNames=getGeneName(gtfFrame)
  for(i in seq(nsize)){
    if(toString(geneNames[i]) %in% namesVec){
      subsetArr[i]=TRUE
    }
  }
  print(sum(subsetArr))
  newDF=gtfFrame[subsetArr,]
  if(length(savename)>0){
    write.table(newDF,file=savename,quote=FALSE, sep='\t', col.names = FALSE,row.names=FALSE)
  }
  return(newDF)
}

#Code to use for deseq style object which has already had rows selected
#Called by subsetGenesDeseq after it chooses genes to retain
#If savename is provided, the resulting gtf will be saved with that name to the current working directory

subsetGeneName = function(gtfFrame,namesVec,savename=""){
  return(subsetGeneNameHelper(gtfFrame,row.names(namesVec),savename=savename))
}

#Function designed to read a file and use the names within to subset a gtf file
#Inputs
#Requires a gtf file, which is in the 10 column vectorbase format
#Requires a file name, where the file is a unix format plaintext file with one gene on each line
#Example below, ignoring #:
#AAEL009762
#AAEL009763
#AAEL009764

#If savename is provided, the resulting gtf will be saved with that name to the current working directory

subsetGeneList = function(gtfFrame,namesFile,savename=""){
  namesVec=readLines(namesFile)
  return(subsetGeneNameHelper(gtfFrame,namesVec,savename=savename))
}

#Function designed to subset a DESeq object and 
#Inputs
#Requires a gtf file, which is in the 10 column vectorbase format
#Requires a DESeq object (or equivalent dataframe with row names that correspond to genes)
#By default, only returns gtf rows with genes that have an adjusted p value of less than 0.05
#This value can be adjusted
#Can also or additionally subset based on the absolute value of log2fc. To toggle an option off, set it negative
#IF ROWS HAVE ALREADY BEEN SELECTED, USE THE subsetGeneName FUNCTION INSTEAD
subsetGenesDeseq = function(gtfFile,deseqOb,pval=0.05,log2fc=-1,savename=""){
  gtfFrame = read_tsv(gtfFile,col_names = FALSE)
  if(pval > 0){
    deseqOb=deseqOb[deseqOb$padj<pval,]
  }
  if(log2fc > 0){
    deseqOb=deseqOb[abs(deseqOb$log2FoldChange)<log2fc,]
  }
  newDF=subsetGeneName(gtfFile,deseqOb,savename=savename)
  return(newDF)
}

#Examples

#Get the gene names in a gtf file in AAEL#### Format
testGN=getGeneName(wholeGTF)

#Save GTF with only the names found in "geneNames.ex" to "ExampleSubset.gtf"
subsetGeneList(wholeGTF,"geneNames.ex",savename="ExampleSubset.gtf")

#Return GTF with only the genes significantly expressed in rna-seq at a p<0.1 level
exDEseqOb = readRDS("RNAseqResults.RDS")
myGTF=subsetGenesDeseq(wholeGTF,exDEseqOb,pval=0.1)

#Return gtf with only 3 genes manually chosen
subsetGeneNameHelper(wholeGTF,c("AAEL009762","AAEL009762","AAEL009762"))

test=read.table("RVFVvBF7dresults_p10_min10-pilot.ex")

wholeGTF=read_tsv("Aedes68GenesExtra.gtf",col_names = FALSE)
wholeGTF=read_tsv("AedesGenes.bed",col_names = FALSE)
wholeGTF=read_tsv("Aedes68.promoters.bed",col_names = FALSE)

subsetGeneList(wholeGTF,"RVFVvBF7dresults_p10_min10-pilot.ex",savename="ExampleSubset2.bed")

subsetGeneNameHelper(wholeGTF,c("AAEL009762","AAEL009762","AAEL009762"),savename = "Test2.bed")

getGeneName(wholeGTF)

substring(wholeGTF[1,10],10,19)

#A function that reads a 10 column gff/gtf file and returns a modified version changed to instead provide regions
#a given number of basepairs upstream of the start site
#First argument is a 10 column gtf file 
#Ex: AedesGenes.gtf
#Second argument is the length

#bedFormat is an optional argument that should be set to "TRUE" if using a bed-style file
#with genomic coordinates in the second and third columns
#It is false by default and the function expects GTF/GFF style files if set to "FALSE"

#geneBody is an optional argument regarding whether to include the gene body as well as the promoter
#By default, it is true and the genomic ranges will have the promoter region added and gene body retained
#If set to "FALSE", only the promoter regions will be returned

#Does not return the new dataframe, but rather writes it to a new file with the same name
#Other than "Prom####" before the suffix, where #### is the chosen length of promoter region

makePromFiles = function(bedFile, promLengths,save=TRUE,bedFormat=FALSE,genebody=TRUE,bidirectional=FALSE){
  
  orientCol=-1
  startCol=-1
  endCol=-1
  #nameCol=-1
  #chromCol=-1
  
  
  if(bedFormat){
    orientCol=6
    startCol=2
    endCol=3
    #chromCol=1
    #nameCol=10
  }
  else{
    #GTF/GFF
    orientCol=7
    startCol=4
    endCol=5
    #chromCol=1
    #nameCol=9
  }
  
  
  wholeBed=read_tsv(bedFile,col_names = FALSE)
  fileroot=substr(bedFile,1,nchar(bedFile)-4)
  geneStart = wholeBed[[startCol]]
  geneEnd = wholeBed[[endCol]]
  count = 1
  for(promLen in promLengths){
    newBed = wholeBed
    for(i in seq(length(geneStart))){
      if(newBed[i,orientCol] == "+"){
        if(genebody==FALSE){
          newBed[i,endCol] = geneStart[i] 
        }
        if(geneStart[i] < promLen){
          newBed[i,startCol] = 0
        }
        else{
          newBed[i,startCol] = geneStart[i] - promLen
        }
        #newBed[i,4]
        if(bidirectional){
          newBed[i,endCol]=newBed[i,endCol]+promLen
        }
      }
      else{
        #
        if(genebody==FALSE){
          newBed[i,startCol] = geneEnd[i]
        }
        #newBed[i,startCol] = geneEnd[i]
        newBed[i,endCol] = geneEnd[i]+promLen
        
        if(bidirectional){
          if(newBed[i,startCol]<promLen){
            newBed[i,startCol]=0
          }
          else{
            newBed[i,startCol]=newBed[i,startCol]-promLen
          }
        }
      }
    }
    biString=""
    if(bidirectional){
      biString="Bi"
    }
    if(genebody){
      outFile = paste(substr(bedFile,1,nchar(bedFile)-4),"GeneANDProm",toString(promLen),biString,substr(bedFile,nchar(bedFile)-3,nchar(bedFile)),sep="")
    }
    else{
      outFile = paste(substr(bedFile,1,nchar(bedFile)-4),"Prom",toString(promLen),biString,substr(bedFile,nchar(bedFile)-3,nchar(bedFile)),sep="") 
    }
    write.table(newBed,file=outFile,quote=FALSE, sep='\t', col.names = FALSE,row.names=FALSE)
    count = count+1
  }
}

#Example Usage

#Create single file 
#Note that even with single file, must use the vector syntax of C(###)
makePromFiles("AedesGenes.gtf",c(500))

#Create multiple files at once
makePromFiles("AedesGenes.gtf",c(750,1500,2000))

#Create file with only promoter instead of promoter and gene body
makePromFiles("AedesGenes.gtf",c(500),genebody=FALSE)

#Create file with bed file instead of a gtf
makePromFiles("AedesGenes.bed",c(500),bedFormat = TRUE)

#Create file that hasa distances on both sides of the TSS but excludes the rest of the gene body
makePromFiles("Aedes68Genes.bed",c(500,2000,5000,10000),bedFormat = TRUE,bidirectional=TRUE,genebody = FALSE)
