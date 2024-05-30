#https://lashlock.github.io/compbio/R_presentation.html
#https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
library(DESeq2)
library(ggplot2)

# directory is where your htseq-count output files are located.
directory <- "/Users/pegaheizad/Desktop/rosenberg-lab/rna-seq-1/"

# samplesFiles is a variable which points to your htseq-count output files,
# sampleFiles <- grep("RNA",list.files(directory),value=TRUE)
sampleFiles <- list.files(directory)

# One for one for your sample type
condition <- c('BF','BF','RVFV','RVFV', 'RVFV')


sampleTable <- data.frame(sampleName = c('BF-1','BF-2','RVFV-1','RVFV-2', 'RVFV-3'),
                          fileName = sampleFiles,
                          condition = condition)
library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
# dds
dds <- DESeq(ddsHTSeq)
res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table
# info on what tests and variables were used to generate results:
mcols(res)$description
# summary of results
summary(res) 
# order results from smallest to largest adjusted p-value
htres <- htres[order(htres$padj),]
head(htres)
#saveRDS(htres,file="hisatDESeqResults.RDS")
# significant data with p value < 0.05
resSig <- subset(res, padj < 0.01)


## the following code blocks are for graphical representations of data
# log2 fold change MA plot
plotMA(res, ylim=c(-2,2))
# shrinked log fold change
# remove the noise associated with log2 fold changes
resLFC <- lfcShrink(dds, coef="condition_RVFV_vs_BF", type="apeglm")
plotMA(resLFC, ylim=c(-2,2))

# plotcount
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
# get which samples are which points
counts(dds)[which.min(res$padj), ]
# do the same plotcount with ggplot (optional in case we want to work with ggplot)
text(x = x_coords, y = y_coords, labels = dds$colnames, pos = 3)
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#effect of transformation on variance 
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

# heat map (we might not need this visual representation if its not helpful)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition")])
rownames(df) <- colnames(dds)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)


#heat map of sample-sample distance
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


## PC plot
plotPCA(vsd, intgroup="condition")
## same thing using ggplot
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, label=rownames(pcaData))) +
  geom_point(size=3) +
  geom_text(size=3, vjust=-0.5, hjust=0.5) +  # Add labels to points with some adjustments
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed()
  