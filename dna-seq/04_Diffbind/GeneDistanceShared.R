library(readr)
library("tools")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(readr.show_col_types = FALSE)


#'gtfFile' is a 9 column GTF file denoting the genes of interest. It can be the entire Aedes genome,
#a subset of GOIs (such as those differentially expressed in RNA-seq), or specific regions (such as 3' UTRs).
#a selection of these files are provided on the CUT&RUN github in dnaseq->04_Diffbind->

#If using the previous 10 column bed-style files, call the function with bedFormat=TRUE (see below)

#'diffPeakFile' is a file containing information about peaks, including chromosome
#location on chromosome, p-value, log2fold change, etc. 
#This is typically generated from Diffbind
#Peak files in bed format are now also supported!


#'inGene' and 'nearEnd' are both parameters that change the set of genes returned
#If 'inGene' is True, rows returned will include peaks inside of a gene body even if they are not within the given
#distance of the end of the gene
#If 'nearEnd' is True, rows returned will include peaks within the given distance of the end of the gene even if 
#they do not overlap with the gene body

#bedFormat is an optional argument that should be set to "TRUE" if using a bed-style file
#with genomic coordinates in the second and third columns
#It is false by default and the function expects GTF/GFF style files if set to "FALSE"

#topN is an optional parameter to only accept the top N peaks and disregard matches for any peak lower in the file
#By default, this behavior is disabled and it is recommend 

#'saveResults' allows the user to save the results as both a CSV and RDS

#sampName is an optional parameter that will prefix the names of any files saved
#Take care as if the name is the same as a previous file THE OLD FILE WILL BE OVERWRITTEN

#'verbose' is an optional parameter that will cause the function to print significantly more text while running
#It is not recommended to be used unless debugging the code


#Outputs Column description (column #)
#Chromosome (1)
#The chromosome the gene is one
#Gene Start (2)
#The 


#Gene Orientation (4)
#+ indicates the gene is on the positive sense strand, - is the negative sense strand
#Peak Number (5)
#This is the number of the row the peak is in on the input file (generally diffbind output)
#5 Prime Distance (6)
#This is the minimum distance from one
#3 Prime Distance
#Peak Inside Gene (8)
#This is a boolean value that indicates whether any part of the peak overlaps with the gene body
#Gene ID (9)
#This is the modern ID used to refer to the gene without ambiguity (i.e. AAEL001234 in Aedes aegypti or gene1234 in Culex tarsalis)
#Gene Name (10)
#This is the gene's "common name", if it possesses one. If it does not, this column will have the gene ID.
#Gene Description (11)
#A full name more complete deiscription of what the gene does, if available. 
#EBI Biotype (12)
#Whether the feature is a protein coding protein, a lncRNA, or another category of region
#This is intended to allow flexibility to incorporate other types of data like exons if desired
#

getGeneDistances=function(gtfFile,diffPeakFile,topN=-1,verbose=FALSE,sampName="ChIP",distance=10000,saveResults=TRUE,inGene=TRUE,nearEnd=TRUE,bedFormat=FALSE){
  
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
    sigResDF=read.csv(diffPeakFile)#,quote="")
    dnaStart=sigResDF$start
    dnaEnd=sigResDF$end
    dnaChrom=sigResDF$chr
    dnaFold=sigResDF$Fold
    
    # Added 4/1/2025
    pvalue=sigResDF$p.value
    fdr=sigResDF$FDR
  }
  else{
    sigResDF=as.data.frame(readr::read_tsv(diffPeakFile,col_names=FALSE),stringsAsFactors=FALSE)
    dnaStart=sigResDF[[2]]
    dnaEnd=sigResDF[[3]]
    dnaChrom=as.character(sigResDF[[1]])
    dnaFold=-sigResDF[[8]]
    pvalue=sigResDF[[9]]
    fdr=sigResDF[[10]]
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
  if(length(rnaSamp)==0){
    print("No matching regions found for the given parameters")
    return()
  }
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
  #print(rnaSamp)
  #print("Here 2")
  reducedDF=data.frame(geneChrom,geneStart,geneEnd,geneOrient)#,rnaTitle)
  colnames(reducedDF) = c("Chromosome","Gene Start", "Gene End","Gene Orientation")#,"Gene Name")
  #print(length(withinDist))
  #print(length(rnaSamp))
  dna5Prime=dna5Prime[withinDist]
  dnaSamp=dnaSamp[withinDist]
  
  over3prime=over3prime[withinDist]
  dna3Prime=dna3Prime[withinDist]
  insideGene=insideGene[withinDist]
  over5prime=over5prime[withinDist]
  reducedDF["Peak Number"]=dnaSamp
  reducedDF["5 Prime Dist"]=dna5Prime
  reducedDF["3 Prime Dist"]=dna3Prime
  reducedDF["Peak Inside Gene"]=insideGene
  
  geneID=rep(0,length(reducedDF$'Peak Start'))
  geneName=rep(0,length(reducedDF$'Peak Start'))
  geneDes=rep(0,length(reducedDF$'Peak Start'))
  geneBiotype=rep(0,length(reducedDF$'Peak Start'))
  
  for(i in seq(length(rnaTitle))){
    #print(rnaTitle[i])
    temp=strsplit(rnaTitle[i],split=';')[[1]]
    #Some do not have a name category, need to handle that
    if(length(temp)==3){
      temp[1]=substring(temp[1],first=4)
      geneID[i]=temp[1]
      #temp[2]=substring(temp[2],first=6)
      geneName[i]="NA"
      temp[2]=substring(temp[2],first=13)
      geneDes[i]=temp[2]
      temp[3]=substring(temp[3],first=13)
      geneBiotype[i]=temp[3]
    }
    else{
      temp[1]=substring(temp[1],first=4)
      geneID[i]=temp[1]
      temp[3]=substring(temp[3],first=13)
      geneDes[i]=temp[3]
      
      
      if(temp[4]=="gene_ebi_biotype=lncRNA"){
        temp[2]=substring(temp[2],first=8)
        geneName[i]=temp[2]
        temp[4]=substring(temp[4],first=18)
        geneBiotype[i]=temp[4]
      }
      else{
        temp[2]=substring(temp[2],first=6)
        geneName[i]=temp[2]
        temp[4]=substring(temp[4],first=13)
        geneBiotype[i]=temp[4]
      }
    }
    #print(temp)
  } 
  
  reducedDF["Gene ID"]=geneID
  reducedDF["Gene Name"]=geneName
  reducedDF["Gene Description"]=geneDes
  reducedDF["EBI Biotype"]=geneBiotype
  

  
  
  foldchange = rep(0,length(dnaSamp))
  for(i in seq(length(dnaSamp))){
    foldchange[i]=dnaFold[dnaSamp[i]]
  }
  
  reducedDF["Fold Change"]=foldchange
  
  # Added 4/1/2025
  peak_pvalues=rep(0,length(dnaSamp))
  peak_fdrs=rep(0,length(dnaSamp))
  
  peakStart=rep(0,length(dnaSamp))
  peakEnd=rep(0,length(dnaSamp))
  for(i in seq(1,length(dnaSamp))){
    peakStart[i]=as.integer(dnaStart[dnaSamp[i]])
    peakEnd[i]=as.integer(dnaEnd[dnaSamp[i]])
    # Added 4/1/2025
    peak_pvalues[i]=as.double(pvalue[dnaSamp[i]])
    peak_fdrs[i]=as.double(fdr[dnaSamp[i]])
  }
  reducedDF["Peak Start"]=peakStart
  reducedDF["Peak End"]=peakEnd
  reducedDF["Over5Prime"]=over5prime
  reducedDF["Over3Prime"]=over3prime
  
  reducedDF["Peak Width"]=reducedDF["Peak End"]-reducedDF["Peak Start"]
  reducedDF["Gene Width"]=reducedDF["Gene End"]-reducedDF["Gene Start"]
  
  peakStartOffset=rep(0,length(reducedDF$'Peak Start'))
  peakEndOffset=rep(0,length(reducedDF$'Peak Start'))

  for(i in seq(length(reducedDF$'Peak Start'))){
    if(reducedDF$'Gene Orientation'[i]=="+"){
      peakStartOffset[i]=reducedDF$'Peak Start'[i]-reducedDF$'Gene Start'[i]
      peakEndOffset[i]=reducedDF$'Peak End'[i]-reducedDF$'Gene Start'[i]
    }
    if(reducedDF$'Gene Orientation'[i]=="-"){
      peakStartOffset[i]=reducedDF$'Gene Start'[i]-reducedDF$'Peak Start'[i]
      peakEndOffset[i]=reducedDF$'Gene Start'[i]-reducedDF$'Peak End'[i]
    }
  }

  reducedDF["Peak Start Offset"]=peakStartOffset
  reducedDF["Peak End Offset"]=peakEndOffset
    
  # Added 4/1/2025
  reducedDF["p-value"]=peak_pvalues
  reducedDF["FDR"]=peak_fdrs
  
  if(topN != -1){
    reducedDF=reducedDF[reducedDF["Peak Number"]<=topN,] 
  }
  if(saveResults){
    saveRDS(reducedDF,paste(sampName,"_DistanceDF.RDS"))
    write.csv(reducedDF,paste(sampName,"_DistanceDF.csv"))
  }
  return(reducedDF)
}


ex=getGeneDistances("Aedes68Genes.bed","PilotSEACR_Results_sig.csv",distance = 5000,inGene=FALSE,bedFormat = TRUE)

#Extracting regions only close to 5 prime end
nrow(ex[ex$`5 Prime Dist`<5000,])

ex=ex[ex$`5 Prime Dist`<5000,]

write.csv(ex,"PilotMergedPeaksWithin5KPromoter.csv")

#Look for matches between the first 10 peaks in a file and genes within 10000 bp of those peaks
ExampleObj = getGeneDistances("AedesGenes.bed","NewPilotPromoterPlusBodyMACS2_Ac_Sig.csv",topN=10,distance = 10000)

#Look for matches between all peaks in a set of diffbind output, using an old bed file
ExampleObj = getGeneDistances("AedesGenes.bed","NewPilotPromoterPlusBodyMACS2_Ac_Sig.csv",bedFormat = TRUE,distance = 10000,inGene=FALSE)

testFrame=read.csv("BF_Ac_D7vsRVFV_Ac_D7_SEACRPilot.csv")
nrow(testFrame)
write.csv(head(testFrame,n=1000),"AssignTests.csv")


#BF_H3K27ac_d7_vs_SF_H3K27ac-DiffBind-out.csv
#AedesGenes_CR.gtf

ex=getGeneDistances("AedesGenes_CR.gtf","Top1000BF_Ac_D7_Diffbind.csv",distance = 2000,inGene=FALSE,bedFormat = TRUE)


ex=getGeneDistances("Aedes68Genes.bed","AssignTests.csv",distance = 5000,inGene=FALSE,bedFormat = TRUE)


write.csv(ex,"ExampleChipDistances.csv")

myPilotResults=getGeneDistances("Aedes68Genes.bed","PilotNew_Results.csv",distance=2000,inGene=TRUE,bedFormat=TRUE)

nearend=getGeneDistances("Aedes68Genes.bed","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\BF_Ago2_D1_Rep2_Control_peaks.narrowPeak",distance = 2000,inGene=FALSE,bedFormat = TRUE)

saveRDS(nearend,"BF_Ago2_D1_Rep2_NearEnd.RDS")

inbody=getGeneDistances("Aedes68Genes.bed","C:\\Users\\hunte\\Desktop\\AltChip\\Peaks\\NewPilotDefault\\BF_Ago2_D1_Rep2_Control_peaks.narrowPeak",distance = 2000,nearEnd = FALSE,inGene=TRUE,bedFormat = TRUE)

saveRDS(inbody,"BF_Ago2_D1_Rep2_InBody.RDS")

notend=inbody[inbody$`5 Prime Dist`>2000,]
notend=notend[notend$`3 Prim Dist`>2000,]#Add e to Prime if reran


unique(notend$`Gene ID`)

unique(nearend$`Gene ID`)

intersect(notend$`Gene ID`,nearend$`Gene ID`)

#Examples using peaks from the sample file
allSamps=read.csv("CUT_RUN_Meta_File_MACS2_0.05_keepdup.csv")
PilotAcSamps=allSamps[allSamps$Factor=="H3K27Ac",]
PilotAcSamps=PilotAcSamps[PilotAcSamps$RiftExperiment==TRUE,]

pilotMergedRVFVAcPeaks=PilotAcSamps[PilotAcSamps$Condition=="RVFV",]$Peaks[1]
pilotMergedBFAcPeaks=PilotAcSamps[PilotAcSamps$Condition=="BF",]$Peaks[1]
rvfvAcDistant=getGeneDistances("Aedes68Genes.bed",pilotMergedRVFVAcPeaks,distance = 200000,nearEnd = TRUE,inGene=FALSE,bedFormat = TRUE)

bfAcDistant=getGeneDistances("Aedes68Genes.bed",pilotMergedBFAcPeaks,distance = 100000,nearEnd = TRUE,inGene=FALSE,bedFormat = TRUE)


testcsv=read.csv("BF_H3K27ac_d7_vs_SF_H3K27ac-DiffBind-out-NEW.csv")

dnaStart=testcsv$start
dnaEnd=testcsv$end
dnaChrom=testcsv$chr
dnaFold=testcsv$Fold

# Added 4/1/2025
pvalue=testcsv$p.value
fdr=testcsv$FDR

sum(testcsv$Conc==0)

sum(dnaFold==0)
sum(is.na(dnaFold))
sum(is.na(dnaFold))
sum(is.na(dnaFold))
sum(is.na(dnaFold))
sum(is.na(pvalue))
sum(is.na(fdr))

table=read.table("AedesGenes.gtf",sep="\t")

tabletigs=table$V1

sum(!(dnaChrom %in% tabletigs))

dnaChrom[!(dnaChrom %in% tabletigs)]

unique(testcsv$)
unique(table$V8)
sum(!(tabletigs %in% dnaChrom))

#9566 problem for some reason