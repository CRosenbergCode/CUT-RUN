

#Applying greylists

customGreyList = function(bam,size=30000,p=0.99,gap=16384){
  cols = c("Chromosome","Organism","Version","Length","Type")
  setwd("C:\\Users\\hunte\\Desktop\\AltChip")
  AedesFrame=as.data.frame(readr::read_delim("AedesChromes.txt",delim=" | ",col_names = cols))
  
  
  
  AedesGenomeGR <- GRanges(
    seqnames = AedesFrame$Chromosome,
    IRanges(1, end=AedesFrame$Length))
  #NEED TO CHANGE
  setwd("C:\\Users\\hunte\\Desktop\\AltChip\\bams")
  AedesGrey <- new("GreyList",genomeRegions=AedesGenomeGR)
  AedesGrey = setRegions(AedesGrey,AedesGenomeGR)
  inputReads = countReads(AedesGrey,bam)
  
  inputReads@counts
  
  thresholdedReads = calcThreshold(inputReads,sampleSize = size,p=p)
  gl = makeGreyList(thresholdedReads,maxGap=gap)
  return(gl)
}


diffbindGreylist = function(samples,group1,group2,greylist){
  setwd("C:\\Users\\hunte\\Desktop\\AltChip")
  
  samps <- read.csv(samples)
  greyDBA = dba(sampleSheet=samps)
  greyDBA = dba.blacklist(greyDBA,blacklist=greylist@regions,greylist=greylist@regions)
  
  greyDBA = dba.count(greyDBA,summits=FALSE)
  
  #consPeaks=dba.peakset(consensus_counts)
  greyDBA <- dba.normalize(greyDBA)
  greyDBA <- dba.contrast(greyDBA,group1=group1, group2=group2,minMembers = 2)
  
  greyDBA <- dba.analyze(greyDBA,method=DBA_ALL_METHODS,bBlacklist = FALSE,bGreylist = FALSE)
  greyDBA <- dba.report(greyDBA, method=DBA_EDGER, contrast = 1, th=1)
  greyDBA = annoGR2DF(greyDBA)
  return(greyDBA)
}


subsetMACS2Peaks = function(inFile,outFile=FALSE,cutoff=4,value="q_value",higher=TRUE){
  col = 0
  cols = c("chromosome","start","end","name","enrichment","strand","signal","p","q","peak")
  
  bedsDF=suppressMessages(as.data.frame(readr::read_tsv(inFile,col_names=cols)))
  
  if(value=="q_value"){
    col = 9
  }
  if(value=="p_value"){
    col = 8
  }
  if(value=="fold_enrichment"){
    col = 7
  }
  if(higher){
    newDF=bedsDF[bedsDF[,col]>cutoff,] 
  }
  else{
    newDF=bedsDF[bedsDF[,col]<cutoff,]
  }
  if(outFile!=FALSE){
    write.table(newDF, file=outFile, quote=FALSE, sep='\t', col.names = NA)
  }
  return(newDF)
}


listSubsets = function(inFiles,cutoffs,value="q_value",higher=TRUE){
  for(file in inFiles){
    print(paste("For",file))
    for(cutoff in cutoffs){
      temp=subsetMACS2Peaks(inFile = file,cutoff=cutoff,value=value,higher=higher)
      print(paste("There are", nrow(temp),"peaks with",value,"greater than",cutoff))
    }
  }
}

