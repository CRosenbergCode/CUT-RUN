#https://www.bioconductor.org/help/course-materials/2016/CSAMA/lab-5-chipseq/Epigenetics.html
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(chipseq)
library(biomaRt)
library(stringr)
library(biomartr)
library(GenomicFeatures)
library(AnnotationDbi)

#Must be run for non-model organisms
options(ucscChromosomeNames=FALSE)
setwd("C:\\Users\\hunte\\Desktop\\AltChip")

#Annotations will be plotted as the 
#Reads and peaks will be plotted in the order provided in the list, with those provided first being plotted above those provided last

#The plotting will be saved to an object, which will then be invoked to use it

#The annotation file should be a GFF3 file that contains not only genes but other features like exons (if plotting transcripts)
#Further file types may be supported in future releases
#This relies upon MACS2 style peak files, which are in BED format with a few additional columns for signal and significance
makeGViz=function(annotationFile="",reads=c(),peaks=c(),geneBodies=TRUE,transcripts=FALSE,chromosome="",orgName="Aedes aegypti",
                  peakNames=c(),alignNames=c()){
  
  customPeakImport=function(file){
    macsPeakCols <- c(signalValue = "numeric", pValue = "numeric",
                              qValue = "numeric", peak = "integer")
    
    macsPeak <- import(file, format = "BED",
                            extraCols = macsPeakCols)
    return(macsPeak)
  }
  
  peakLists=c()
  if(geneBodies){
    
    tempgff=read_gff(annotationFile)
    
    geneNames = unlist(tempgff$attribute)
    
    #Currently only for Aedes
    geneIDS=lapply(geneNames,str_extract,pattern="AAEL[0-9]+")
    geneIDS=unlist(geneIDS)
    geneBodyTrack <- AnnotationTrack(range=annotationFile,group=geneIDS,fontcolor.item="black",just.group = "above",name="Aedes aegypti Genes",cex.title=1)
    peakLists=c(peakLists, geneBodyTrack)
    #append(peakList,geneBodyTrack)
  }
  if(transcripts){
    testtx=makeTxDbFromGFF(annotationFile,format="gff3")
    transcriptTrack <- GeneRegionTrack(
      range = testtx,
      genome = "Aedes aegypti",
      name="Aedes aegypti Transcripts",
      cex.title=1
    )
    peakLists=c(peakLists, transcriptTrack)
    #append(peakLists,transcriptTrack)
  }
  if(length(reads>0)){
    if(length(alignNames)==0){
      alignNames=paste0("Alignment_", 1:length(reads))
    }
    for(i in seq(length(reads))){
      alTrack <- DataTrack(range=reads[i],type="histogram",
                           name=alignNames[i],isPaired = TRUE,cex.title=1)
      peakLists=c(peakLists, alTrack)
      #append(peakLists,transcriptTrack)
    }
  }

  if(length(peaks>0)){
    if(length(peakNames)==0){
      peakNames=paste0("Peaks_", 1:length(peaks))
    }
    for(i in seq(length(peaks))){
      tempPeak= AnnotationTrack(peaks[i],#peaks.rep1, 
                                name=peakNames[i],
                                #chromosome='AaegL5_2',
                                importFunction = customPeakImport,
                                shape='box',fill='blue3',size=20,cex.title=1)
      peakLists=c(peakLists, tempPeak)
      #append(peakList,transcriptTrack)
    }
  }

  return(peakLists)
}


#This function is called using the ouput of the makeGviz function (a vector of tracks)
#A chromosome and start and end point must also be provided
plotAligns=function(tracks,chromosome,start,end,gene="",heights=c(),title="",gff="VectorBase-68_AaegyptiLVP_AGWG.gff"){
  if(nchar(gene)>0){
    testtx=makeTxDbFromGFF(gff,format="gff3")
    myannot=AnnotationDbi::select(testtx, keys = c(gene), columns=c("TXCHROM","TXSTART","TXEND","TXNAME"),keytype="GENEID")
    chromosome=myannot$TXCHROM
    start=myannot$TXSTART
    end=myannot$TXEND
  }
  if(length(heights)==0){
    heights=rep(1,length(tracks))
  }
  if(length(title)>0){
    plotTracks(tracks, sizes=heights,main=title,
               from=start, to=end,transcriptAnnotation="symbol",groupAnnotation = "group",just.group = "above",chromosome = chromosome)#,
    
  }
  else{
    plotTracks(tracks, sizes=heights,
               from=start, to=end,transcriptAnnotation="symbol",groupAnnotation = "group",just.group = "above",chromosome = chromosome)#,
    
  }
  #The track height can be controlled by providing a vector of relative heights to the sizes parameter of the plotTracks() function.
  # plotTracks(tracks, sizes=heights,
  #            from=start, to=end,transcriptAnnotation="symbol",groupAnnotation = "group",just.group = "above",chromosome = chromosome)#,
  # 
  
  #plotTracks(c(bfAcTrack,atrack,peaks1.track),
  #           from=start, to=end,transcriptAnnotation="symbol",showId=TRUE,showFeatureId=TRUE,chromosome = chromosome)#,
  
  
  
  #plotTracks(tracks,from=start, to=end,
  #           transcriptAnnotation="symbol",groupAnnotation = "group",just.group = "above")#,
  
  
}

#Examples Below
allSamps=read.csv("CUT_RUN_Meta_File_MACS2_0.05_keepdup.csv")
pilotAcSamps=allSamps[allSamps$RiftExperiment==TRUE,]
pilotAcSamps=pilotAcSamps[pilotAcSamps$Factor=="H3K27Ac",]
pilotAcSamps=pilotAcSamps[pilotAcSamps$Day ==7,]
pilotAcBFPeaks=pilotAcSamps[3,]$Peaks
pilotAcRVFVPeaks=pilotAcSamps[4,]$Peaks

acTracks=pilotAcSamps[pilotAcSamps$Replicate!=3,]$bamReads

#Get the trac
mytracks=makeGViz(annotationFile="VectorBase-68_AaegyptiLVP_AGWG.gff",peaks=c(pilotAcBFPeaks,pilotAcRVFVPeaks),geneBodies=TRUE,transcripts=TRUE,peakNames=c("RVFV_D7_Ac","BF_D7_Ac"),
                  reads=acTracks,alignNames = c("BF_Ac_D7_1","BF_Ac_D7_2","RVFV_Ac_D7_1.1","RVFV_Ac_D7_1.2","RVFV_Ac_D7_2"))



mergedAcTracks=c("bams/RenamedBams/BF_H3K27Ac_Merged_Aae.bam","bams/RenamedBams/RVFV_H3K27Ac_Merged_Aae.bam")
myMergedTracks=makeGViz(annotationFile="VectorBase-68_AaegyptiLVP_AGWG.gff",peaks=c(pilotAcBFPeaks,pilotAcRVFVPeaks),geneBodies=TRUE,transcripts=TRUE,peakNames=c("BF_D7_Ac","RVFV_D7_Ac"),
                  reads=mergedAcTracks,alignNames = c("BF_Ac_Merged","RVFV_Ac_Merged"))



mytracks=makeGViz(annotationFile="VectorBase-68_AaegyptiLVP_AGWG.gff",peaks=c(pilotAcBFPeaks,pilotAcRVFVPeaks),geneBodies=TRUE,peakNames=c("RVFV_D7_Ac","BF_D7_Ac"),
                  reads=c())

plotAligns(mytracks,chromosome="AaegL5_2",start=4060221, end=4134232)

plotAligns(mytracks,chromosome="AaegL5_2",start=4080221, end=4114232,heights=c(2,2,2,2,2,2,2,1,1))


plotAligns(myMergedTracks,chromosome="AaegL5_2",start=4060221, end=4134232)



mytracks=makeGViz(annotationFile="VectorBase-68_AaegyptiLVP_AGWG.gff",peaks=c(pilotAcBFPeaks,pilotAcRVFVPeaks),geneBodies=TRUE,transcripts=TRUE,peakNames=c("RVFV_Peaks","BF_Peaks"),
                  reads=acTracks,alignNames = c("BF_1","BF_2","RVFV_1.1","RVFV_1.2","RVFV_2"))

#Get the coordinates of genes
testtx=makeTxDbFromGFF("VectorBase-68_AaegyptiLVP_AGWG.gff",format="gff3")
mykeys=c("AAEL009762","AAEL004591","AAEL019495","AAEL021746","AAEL005849","AAEL007050","AAEL023746")
AnnotationDbi::select(testtx, keys = mykeys, columns=c("TXCHROM","TXSTART","TXEND","TXNAME"),keytype="GENEID")

# AAEL009762	cytochrome P450
# AAEL004591	Dbl homology (DH) domain, RHO GUANINE NUCLEOTIDE EXCHANGE FACTOR 
# AAEL019495	Pleckstrin homology domain;PDZ domain
# AAEL021746	TLDc domain
# AAEL005849	synaptic vesicle protein
# AAEL007050	sugar transporter
# AAEL023746	Leucine-rich repeat N-terminal domain;Leucine-rich repeat

#Plot list of genes of interest
plotAligns(myMergedTracks,chromosome="AaegL5_1",start=4080222, end=4134232,title="AAEL009762 cytochrome P450")

plotAligns(myMergedTracks,chromosome="AaegL5_2",start=32872874, end=33496772,title="AAEL004591")

plotAligns(myMergedTracks,chromosome="AaegL5_2",start=418526388, end=419058148,title="AAEL019495")

tempSave=plotAligns(myMergedTracks,chromosome="AaegL5_3",start=124918308, end=125554053,title="AAEL021746")

plotAligns(myMergedTracks,chromosome="AaegL5_3",start=218540137, end=218633317,title="AAEL005849")

plotAligns(myMergedTracks,chromosome="AaegL5_3",start=298721843, end=298756329,title="AAEL007050")

plotAligns(myMergedTracks,chromosome="AaegL5_3",start=342318151, end=342385526,title="AAEL023746")



#Plot using input subtracted bigwig files
bwAcTracks=c("bams/RenamedBams/compare_Merged_BF_Ac12.bw","bams/RenamedBams/compare_Merged_RVFV12_Ac.bw")


bwMergedTracks=makeGViz(annotationFile="VectorBase-68_AaegyptiLVP_AGWG.gff",peaks=c(pilotAcBFPeaks,pilotAcRVFVPeaks),geneBodies=TRUE,transcripts=TRUE,peakNames=c("BF_D7_Ac","RVFV_D7_Ac"),
                        reads=bwAcTracks,alignNames = c("BF_Ac_Merged","RVFV_Ac_Merged"))

plotAligns(bwMergedTracks,chromosome="AaegL5_1",start=4080222, end=4134232,title="AAEL009762 cytochrome P450")

plotAligns(bwMergedTracks,chromosome="AaegL5_2",start=32872874, end=33496772,title="AAEL004591")

plotAligns(bwMergedTracks,chromosome="AaegL5_2",start=418526388, end=419058148,title="AAEL019495")

tempSave=plotAligns(bwMergedTracks,chromosome="AaegL5_3",start=124918308, end=125554053,title="AAEL021746")

plotAligns(bwMergedTracks,chromosome="AaegL5_3",start=218540137, end=218633317,title="AAEL005849")

plotAligns(bwMergedTracks,chromosome="AaegL5_3",start=298721843, end=298756329,title="AAEL007050")

plotAligns(bwMergedTracks,chromosome="AaegL5_3",start=342318151, end=342385526,title="AAEL023746")