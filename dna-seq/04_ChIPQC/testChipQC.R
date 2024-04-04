#https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/05_combine_chipQC_and_metrics.html
library(ChIPQC)
setwd("C:\\Users\\hunte\\Desktop\\ChipData")
#BiocManager::install("ChIPQC")
samples <- read.csv('chipSamples2_28_24.csv')
View(samples)
register(SerialParam())


chipObj <- ChIPQC(samples)#, annotation="hg19") 
ChIPQCreport(chipObj, reportName="ChIP QC report: Mosquito K9", reportFolder="ChIPQCreport")
chipObj

length(samples)
head(samples)
unique(c(rownames(meta), names(controlist)))

samplesIn <- read.csv('chipSamplesWInput.csv')
head(samplesIn)
chipObjIn <- ChIPQC(samplesIn)#, annotation="hg19") 
head(chipObjIn)
reportIn = ChIPQCreport(chipObjIn, reportName="ChIP QC report: CUTRUN",facetBy = "Condition",colourBy="Factor", reportFolder="ChIPQCreport")
reportIn

test = data(tamoxifen_QC)

test = ChIPQCreport(tamoxifen,facetBy="Tissue",colourBy="Condition")

