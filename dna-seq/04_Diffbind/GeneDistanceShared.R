library(readr)
setwd("C:\\Users\\hunte\\Desktop\\ChipData")
options(readr.show_col_types = FALSE)

getGeneDistances=function(bedFile,diffPeakFile,topN=FALSE,verbose=FALSE,sampName="ChIP",distance=1000,saveResults=TRUE,inGene=TRUE,nearEnd=TRUE){
  
  rnaDF=as.data.frame(readr::read_tsv(bedFile,col_names=FALSE),stringsAsFactors=FALSE,show_co)
  rnaNames=as.list(rnaDF[10][,])
  length(rnaNames)

  sigResDF=read.csv("AcvsOthersDiffPeaks.csv")#diffPeakFile)

  rnaStart=as.list(rnaDF[2][,])
  rnaEnd=as.list(rnaDF[3][,])
  rnaChrom=as.list(rnaDF[1][,])
  
  dnaStart=as.list(sigResDF[2][,])
  dnaEnd=as.list(sigResDF[3][,])
  #print(dnaStart)
  #print(dnaEnd)
  
  dnaChrom=as.character(sigResDF[1][,])
  
  withinDist=rep(FALSE,length(rnaChrom)*length(dnaStart))
  dnaDist=rep(0,length(rnaChrom)*length(dnaStart))
  dna5Prime=rep(0,length(rnaChrom)*length(dnaStart))
  dnaSamp=rep(0,length(rnaChrom)*length(dnaStart))
  insideGene=rep(FALSE,length(rnaChrom)*length(dnaStart))
  rnaSamp=rep(0,length(rnaChrom)*length(dnaStart))
  over5prime=rep(FALSE,length(rnaChrom)*length(dnaStart))
  
  rnaLen=length(rnaChrom)
  dnaLen=length(dnaStart)
  for(i in 1:rnaLen){
    orientation=rnaDF[i,6]
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
        if(inGene && ((as.integer(rnaEnd[i]) < as.integer(dnaEnd[j])) && (as.integer(rnaStart[i]) > as.integer(dnaStart[j])) || (as.integer(rnaStart[i]) < as.integer(dnaStart[j])) && (as.integer(rnaEnd[i]) > as.integer(dnaEnd[j])) || (as.integer(rnaEnd[i]) > as.integer(dnaStart[j])) && (as.integer(rnaEnd[i]) < as.integer(dnaEnd[j])) || (as.integer(rnaStart[i]) > as.integer(dnaStart[j])) && (as.integer(rnaStart[i]) < as.integer(dnaEnd[j])))){
          if(verbose){
            print("Inside the gene!")
          }
          
          if((as.integer(rnaEnd[i]) > as.integer(dnaStart[j])) && (as.integer(rnaEnd[i]) < as.integer(dnaEnd[j])) && orientation=='-'){
            over5prime[tempNum]=TRUE
          }
          if((as.integer(rnaStart[i]) > as.integer(dnaStart[j])) && (as.integer(rnaStart[i]) < as.integer(dnaEnd[j])) && orientation=='+'){
            over5prime[tempNum]=TRUE
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
          }
          if(orientation=='-'){
            dna5Prime[tempNum]=min(c(abs(as.integer(rnaEnd[i])-as.integer(dnaEnd[j])),abs(as.integer(rnaEnd[i])-as.integer(dnaStart[j]))))
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
    geneChrom[i]=rnaDF[k,1]
    geneOrient[i]=rnaDF[k,6]
    geneStart[i]=rnaDF[k,2]
    geneEnd[i]=rnaDF[k,3]
    rnaTitle[i]=rnaDF[k,10]
  }

  #print("Here 2")
  reducedDF=data.frame(geneChrom,geneStart,geneEnd,geneOrient,rnaTitle)
  colnames(reducedDF) = c("Chromosome","Gene Start", "Gene End","Gene Orientation","Gene Name")
  #print(length(withinDist))
  #print(length(rnaSamp))
  dna5Prime=dna5Prime[withinDist]
  dnaSamp=dnaSamp[withinDist]
  #print(length(dnaSamp))
  dnaDist=dnaDist[withinDist]
  insideGene=insideGene[withinDist]
  over5prime=over5prime[withinDist]
  reducedDF["5 Prime Dist"]=dna5Prime
  reducedDF["Peak Number"]=dnaSamp
  reducedDF["Gene Dist"]=dnaDist
  reducedDF["Peak Inside Gene"]=insideGene
  #print("Here 4")
  #print(length(dnaSamp))
  peakStart=rep(0,length(dnaSamp))
  peakEnd=rep(0,length(dnaSamp))
  for(i in seq(1,length(dnaSamp))){
    peakStart[i]=as.integer(dnaStart[dnaSamp[i]])
    peakEnd[i]=as.integer(dnaEnd[dnaSamp[i]])
  }
  reducedDF["Peak Start"]=peakStart
  reducedDF["Peak End"]=peakEnd
  reducedDF["Over5Prime"]=over5prime
  #print(reducedDF)
  #print("Here 3")
  reducedDF=reducedDF[reducedDF["Peak Number"]<=topN,]
  if(saveResults){
    saveRDS(reducedDF,paste(sampName,"_DistanceDF.RDS"))
  }
  return(reducedDF)
}

#grep("AAEL[0-9]+",sample)


#test = getGeneDistances("AedesGenes.bed","macs2EDGER_Ac.csv",topN=200,verbose=FALSE)
#getGeneDistances("sig48RNA.bed","seacrEDGER_Ac.csv",verbose=TRUE,distance=10000)
#test = getGeneDistances("AedesGenes.bed","macs2EDGER_Ac.csv",topN=400,verbose=FALSE)
#seacrTest = getGeneDistances("AedesGenes.bed","seacrEDGER_Ac.csv",topN=700,verbose=FALSE,inGene = FALSE,distance=1e4)

SpencerFrame2 = getGeneDistances("AedesGenes.bed","AcvsOthersDiffPeaks.csv",topN=100,verbose=FALSE)

test = SpencerFrame2[,5]

for(i in test){
  #print("Test")
  #print(i)
  print(grep("AAEL[0-9]+",i,value=TRUE))
}

#library(stringr)
#str_match(alice, ".*\\.D([0-9]+)\\.LIS.*")[, 2]