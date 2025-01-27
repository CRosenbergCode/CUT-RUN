library(readr)
library("tools")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(readr.show_col_types = FALSE)

getGeneDistances=function(gtfFile,diffPeakFile,topN=-1,verbose=FALSE,sampName="ChIP",distance=1000,saveResults=TRUE,inGene=TRUE,nearEnd=TRUE,bedFormat=FALSE){
  
  orientCol=-1
  startCol=-1
  endCol=-1
  nameCol=-1
  chromCol=-1
  
  
  if(bedFormat){
    orientCol=6
    startCol=2
    endCol=3
    chromCol=1
    nameCol=10
  }
  else{
    #GTF/GFF
    orientCol=7
    startCol=4
    endCol=5
    chromCol=1
    nameCol=9
  }
  
  rnaDF=as.data.frame(readr::read_tsv(gtfFile,col_names=FALSE),stringsAsFactors=FALSE,show_co)
  rnaNames=as.list(rnaDF[nameCol][,])
  
  diffExt=file_ext(diffPeakFile)

  sigResDF=""
  

  if(diffExt == "csv"){
    sigResDF=read.csv(diffPeakFile)
    dnaStart=sigResDF$start
    dnaEnd=sigResDF$end
    dnaChrom=sigResDF$chr
    dnaFold=-sigResDF$Fold
  }
  else{
    sigResDF=as.data.frame(readr::read_tsv(diffPeakFile,col_names=FALSE),stringsAsFactors=FALSE)
    dnaStart=sigResDF[[2]]
    dnaEnd=sigResDF[[3]]
    dnaChrom=as.character(sigResDF[[1]])
    dnaFold=sigResDF[[8]]
    #dnaStart=sigResDF$start
    #dnaEnd=sigResDF$end
    #dnaChrom=sigResDF$chr
  }
  
  rnaStart=rnaDF[[startCol]]
  rnaEnd=rnaDF[[endCol]]
  rnaChrom=rnaDF[[chromCol]]
  
  #dnaStart=sigResDF$start
  #dnaEnd=sigResDF$end
  #dnaChrom=sigResDF$chr
  
  withinDist=rep(FALSE,length(rnaChrom)*length(dnaStart))
  dnaDist=rep(0,length(rnaChrom)*length(dnaStart))
  dna5Prime=rep(0,length(rnaChrom)*length(dnaStart))
  dnaSamp=rep(0,length(rnaChrom)*length(dnaStart))
  insideGene=rep(FALSE,length(rnaChrom)*length(dnaStart))
  rnaSamp=rep(0,length(rnaChrom)*length(dnaStart))
  over5prime=rep(FALSE,length(rnaChrom)*length(dnaStart))
  
  #New Additions
  over3prime=rep(FALSE,length(rnaChrom)*length(dnaStart))
  dna3Prime=rep(0,length(rnaChrom)*length(dnaStart))
  #geneWidth=rep(0,length(rnaChrom)*length(dnaStart))
  #over5prime=rep(FALSE,length(rnaChrom)*length(dnaStart))
  rnaLen=length(rnaChrom)
  dnaLen=length(dnaStart)
  for(i in 1:rnaLen){
    orientation=rnaDF[i,orientCol]
    for(j in 1:dnaLen){
      tempNum=i*dnaLen+j
      if(dnaChrom[j]==rnaChrom[i]){
        setCondition=FALSE
        helpPrintFunc=function(){
          print(paste("On Chromosome ",dnaChrom[j]))
          print(paste("DNA Peak number",j))
          print(paste("Start of DNA Peak:",toString(dnaStart[j])))
          print(paste("End of DNA Peak:",toString(dnaEnd[j])))
          print(rnaNames[i])
          print(paste("Start of Gene:",toString(rnaStart[i])))
          cat(paste("End of Gene",toString(rnaEnd[i]),"\n\n"))
        }
        dist=Inf
        temp=abs(as.integer(rnaStart[i])-as.integer(dnaStart[j]))
        if(temp < dist){
          dist = temp
        }
        temp=abs(as.integer(rnaStart[i])-as.integer(dnaEnd[j]))
        if(temp < dist){
          dist = temp
        }
        temp=abs(as.integer(rnaEnd[i])-as.integer(dnaStart[j]))
        if(temp < dist){
          dist = temp
        }
        temp=abs(as.integer(rnaEnd[i])-as.integer(dnaEnd[j]))
        if(temp < dist){
          dist = temp
        }
        if((nearEnd) && (dist<distance)){#1e3){
          if(verbose){
            print("RNA")
            print(rnaNames[i])
            #print(i)
            print("Distance is")
            print(dist)
          }
          setCondition=TRUE

        }
        #All four possible scenarios for overlapping with gene body
        if(inGene && ((as.integer(rnaEnd[i]) < as.integer(dnaEnd[j])) && (as.integer(rnaStart[i]) > as.integer(dnaStart[j])) || (as.integer(rnaStart[i]) < as.integer(dnaStart[j])) && (as.integer(rnaEnd[i]) > as.integer(dnaEnd[j])) || (as.integer(rnaEnd[i]) > as.integer(dnaStart[j])) && (as.integer(rnaEnd[i]) < as.integer(dnaEnd[j])) || (as.integer(rnaStart[i]) > as.integer(dnaStart[j])) && (as.integer(rnaStart[i]) < as.integer(dnaEnd[j])))){
          if(verbose){
            print("Inside the gene!")
          }
          #Condition for when the peak overlap 5' end of gene on '-' strand
          if((as.integer(rnaEnd[i]) > as.integer(dnaStart[j])) && (as.integer(rnaEnd[i]) < as.integer(dnaEnd[j])) && orientation=='-'){
            over5prime[tempNum]=TRUE
          }
          #Condition for when the peak overlap the 5' end of the gene on the '+' strand
          if((as.integer(rnaStart[i]) > as.integer(dnaStart[j])) && (as.integer(rnaStart[i]) < as.integer(dnaEnd[j])) && orientation=='+'){
            over5prime[tempNum]=TRUE
          }
          
          #Condition for when the peak overlap the 3' end of the gene on the '-' strand
          if((as.integer(rnaStart[i]) > as.integer(dnaStart[j])) && (as.integer(rnaStart[i]) < as.integer(dnaEnd[j])) && orientation=='-'){
            over3prime[tempNum]=TRUE
          }
          
          
          #Condition for when the peak overlap the 3' end of the gene on the '+' strand
          if((as.integer(rnaEnd[i]) > as.integer(dnaStart[j])) && (as.integer(rnaEnd[i]) < as.integer(dnaEnd[j])) && orientation=='+'){
            over3prime[tempNum]=TRUE
          }
          
          insideGene[tempNum]=TRUE
          setCondition=TRUE
          #withinDist[i]=TRUE
        }
        if(setCondition==TRUE){
          if(verbose){
            helpPrintFunc()
          }
          withinDist[tempNum]=TRUE
          if(orientation=='+'){
            dna5Prime[tempNum]=min(c(abs(as.integer(rnaStart[i])-as.integer(dnaEnd[j])),abs(as.integer(rnaStart[i])-as.integer(dnaStart[j]))))
            #New
            dna3Prime[tempNum]=min(c(abs(as.integer(rnaEnd[i])-as.integer(dnaEnd[j])),abs(as.integer(rnaEnd[i])-as.integer(dnaStart[j]))))
          }
          if(orientation=='-'){
            dna5Prime[tempNum]=min(c(abs(as.integer(rnaEnd[i])-as.integer(dnaEnd[j])),abs(as.integer(rnaEnd[i])-as.integer(dnaStart[j]))))
            #New
            dna3Prime[tempNum]=min(c(abs(as.integer(rnaStart[i])-as.integer(dnaEnd[j])),abs(as.integer(rnaStart[i])-as.integer(dnaStart[j]))))
          }
          dnaSamp[tempNum]=j
          dnaDist[tempNum]=dist
          rnaSamp[tempNum]=i
        }
      }
    }
  }
  rnaSamp=rnaSamp[withinDist]
  geneChrom=rep(0,length(rnaSamp))
  geneOrient=rep(0,length(rnaSamp))
  geneStart=rep(0,length(rnaSamp))
  geneEnd=rep(0,length(rnaSamp))
  rnaTitle=rep(0,length(rnaSamp))
  for(i in seq(1:length(rnaSamp))){
    k=rnaSamp[i]
    #print(i)
    geneChrom[i]=rnaDF[k,chromCol]
    geneOrient[i]=rnaDF[k,orientCol]
    geneStart[i]=rnaDF[k,startCol]
    geneEnd[i]=rnaDF[k,endCol]
    rnaTitle[i]=rnaDF[k,nameCol]
  }

  #print("Here 2")
  reducedDF=data.frame(geneChrom,geneStart,geneEnd,geneOrient,rnaTitle)
  colnames(reducedDF) = c("Chromosome","Gene Start", "Gene End","Gene Orientation","Gene Name")
  #print(length(withinDist))
  #print(length(rnaSamp))
  dna5Prime=dna5Prime[withinDist]
  dnaSamp=dnaSamp[withinDist]
  #print(length(dnaSamp))
  #dnaDist=dnaDist[withinDist]
  over3prime=over3prime[withinDist]
  dna3Prime=dna3Prime[withinDist]
  insideGene=insideGene[withinDist]
  over5prime=over5prime[withinDist]
  reducedDF["Peak Number"]=dnaSamp
  reducedDF["5 Prime Dist"]=dna5Prime
  reducedDF["3 Prim Dist"]=dna3Prime
  reducedDF["Peak Inside Gene"]=insideGene

  
  foldchange = rep(0,length(dnaSamp))
  for(i in seq(length(dnaSamp))){
    foldchange[i]=dnaFold[dnaSamp[i]]
  }
  
  reducedDF["Fold Change"]=foldchange
  
  peakStart=rep(0,length(dnaSamp))
  peakEnd=rep(0,length(dnaSamp))
  for(i in seq(1,length(dnaSamp))){
    peakStart[i]=as.integer(dnaStart[dnaSamp[i]])
    peakEnd[i]=as.integer(dnaEnd[dnaSamp[i]])
  }
  reducedDF["Peak Start"]=peakStart
  reducedDF["Peak End"]=peakEnd
  reducedDF["Over5Prime"]=over5prime
  reducedDF["Over3Prime"]=over3prime
  
  reducedDF["Peak Width"]=reducedDF["Peak End"]-reducedDF["Peak Start"]
  reducedDF["Gene Width"]=reducedDF["Gene End"]-reducedDF["Gene Start"]
  
  #print(reducedDF)
  #print("Here 3")
  if(topN != -1){
    reducedDF=reducedDF[reducedDF["Peak Number"]<=topN,] 
  }
  if(saveResults){
    saveRDS(reducedDF,paste(sampName,"_DistanceDF.RDS"))
    write.table(reducedDF,paste(sampName,"_DistanceDF.csv"))
  }
  return(reducedDF)
}


ExampleObj = getGeneDistances("AedesGenes.bed","Me_RVFV_minus_BF.broadPeak",topN=10000,verbose=FALSE)

ExampleObj = getGeneDistances("AedesGenes.bed","AcvsOthersDiffPeaks.csv",verbose=FALSE)

ExampleObj = getGeneDistances("Aedes68Genes.bed","NewPilotPromoterPlusBodyMACS2_Ac_Sig.csv",verbose=FALSE,distance = 10000,bedFormat = TRUE)

setwd("C:\\Users\\hunte\\Desktop\\ChipData")
setwd("C:\\Users\\hunte\\Desktop\\AltChip")

#write.csv(ExampleObj,"PilotPromoterPlusBodyMACS2_Ac_GeneDistances.csv")

ExampleObj = getGeneDistances("Aedes68GenesExtra.gtf","NewPilotPromoterPlusBodyMACS2_Ac_Sig.csv",verbose=FALSE,distance = 10000)


#'gtfFile' is a 9 column GTF file denoting the genes of interest. It can be the entire Aedes genome,
#or a subset of GOIs (such as those differentially expressed in RNA-seq).

#'diffPeakFile' is a file containing information about peaks, including chromosome
#location on chromosome, p-value, log2fold change, etc. 
#This is typically generated from Diffbind
#Peak files in bed format are now also supported!

#'samples1' and 'samples2' are indices which indicate the samples to compare. 
#The numbers are equivalent to rows in the csv or dataframe, excluding headers
#As a reminder, R is 1 indexed so the first object will be 1

#'plotting' simply determines whether or not plots will be displayed
#It is turned off by default

#'saving' allows the user to save intermediate files such as 
#Take care as if the name is the same as a previous file THE OLD FILE WILL BE OVERWRITTEN
#If the file summary file is to be saved, users should do this themselves.

#