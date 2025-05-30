Genes2000=read.table("Genes2000.txt")

Genes1000=read.table("Genes1000.txt")

library(org.Aaegypti.eg.db)

ego2000 <- enrichGO(gene = unlist(Genes2000), 
                    keyType = "GID", 
                    OrgDb = org.Aaegypti.eg.db, 
                    ont = "CC", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.9,
                    readable = TRUE,minGSSize = 10,maxGSSize = 100)
dotplot(ego2000)


ego1000 <- enrichGO(gene = unlist(Genes1000), 
                    keyType = "GID", 
                    OrgDb = org.Aaegypti.eg.db, 
                    ont = "All", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.9,
                    readable = TRUE,minGSSize = 20,maxGSSize = 500)
dotplot(ego1000)


peaks_RSF1=read.table("peaks_RSF1_RVFV.txt")

mappingPeaks=read.csv("PilotAcRVFV50K_200K_Distal.csv")



#mappingPeaks$V4 %in% unlist(peaks_RSF1)
#unlist(peaks_RSF1) %in% mappingPeaks$V4
#testPeaks=mappingPeaks[which(mappingPeaks$V4==peaks_RSF1),]
testPeaks=mappingPeaks[mappingPeaks$V4 %in% unlist(peaks_RSF1),]

#mappingPeaks$V4 %in% unlist(peaks_RSF1)
testIDs=unique(testPeaks$geneId)

#peaks_RSF1_RVFV.txt


egoRSF1 <- enrichGO(gene = testIDs, 
                    keyType = "GID", 
                    OrgDb = org.Aaegypti.eg.db, 
                    ont = "All", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.9,
                    readable = TRUE)#,minGSSize = 20,maxGSSize = 500)
dotplot(egoRSF1)





egoRSF1 <- enrichGO(gene = testIDs, 
                    keyType = "GID", 
                    OrgDb = org.Aaegypti.eg.db, 
                    ont = "MF", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.9,
                    readable = TRUE)#,minGSSize = 20,maxGSSize = 500)
dotplot(egoRSF1)


egoRSF1MF <- enrichGO(gene = testIDs, 
                      keyType = "GID", 
                      OrgDb = org.Aaegypti.eg.db, 
                      ont = "MF", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.1,
                      readable = TRUE)#,minGSSize = 20,maxGSSize = 500)
dotplot(egoRSF1MF)

edox2 <- pairwise_termsim(egoRSF1MF)
p1 <- treeplot(edox2,nCluster = 4)
p1

egoRSF1CC <- enrichGO(gene = testIDs, 
                      keyType = "GID", 
                      OrgDb = org.Aaegypti.eg.db, 
                      ont = "CC", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.9,
                      readable = TRUE)#,minGSSize = 20,maxGSSize = 500)
dotplot(egoRSF1CC)



egoRSF1BP <- enrichGO(gene = testIDs, 
                      keyType = "GID", 
                      OrgDb = org.Aaegypti.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.1,
                      readable = TRUE)#,minGSSize = 20,maxGSSize = 500)
dotplot(egoRSF1BP)


egoRSF1BPSelect <- enrichGO(gene = testIDs, 
                            keyType = "GID", 
                            OrgDb = org.Aaegypti.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            pvalueCutoff = 0.9,
                            readable = TRUE,minGSSize = 10,maxGSSize = 100)
dotplot(egoRSF1BPSelect)


myresult=egoRSF1BPSelect@result






peaks_ADF1=read.table("ADF1Genes.txt")



egoADF1 <- enrichGO(gene = unlist(peaks_ADF1), 
                    keyType = "GID", 
                    OrgDb = org.Aaegypti.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 1,
                    readable = TRUE,minGSSize = 10,maxGSSize = 40)
dotplot(egoADF1)



peaks_CDC5=read.table("CDC5Genes.txt")

egoCDC5 <- enrichGO(gene = unlist(peaks_CDC5), 
                    keyType = "GID", 
                    OrgDb = org.Aaegypti.eg.db, 
                    ont = "MF", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 1,
                    readable = TRUE,minGSSize = 10,maxGSSize = 100)
dotplot(egoCDC5)



peaks_Hoxc9=read.table("peaks_Hoxc9.txt")

mappingPeaks=read.csv("PilotAcRVFV50K_200K_Distal.csv")


#mappingPeaks$V4 %in% unlist(peaks_RSF1)
#unlist(peaks_RSF1) %in% mappingPeaks$V4
#testPeaks=mappingPeaks[which(mappingPeaks$V4==peaks_RSF1),]


egoMap <- enrichGO(gene = mappingPeaks$geneId, 
                   keyType = "GID", 
                   OrgDb = org.Aaegypti.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.9,
                   readable = TRUE)#,minGSSize = 20,maxGSSize = 500)
dotplot(egoMap)
edox2 <- pairwise_termsim(egoMap)
p1 <- treeplot(edox2)
p1

#mappingPeaks$V4 %in% unlist(peaks_RSF1)
testIDs=unique(testPeaks$geneId)

#peaks_RSF1_RVFV.txt

egoHoxc9 <- enrichGO(gene = testIDs, 
                     keyType = "GID", 
                     OrgDb = org.Aaegypti.eg.db, 
                     ont = "All", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 0.9,
                     readable = TRUE)#,minGSSize = 20,maxGSSize = 500)
dotplot(egoHoxc9)



egoHoxc9 <- enrichGO(gene = testIDs, 
                     keyType = "GID", 
                     OrgDb = org.Aaegypti.eg.db, 
                     ont = "All", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 0.9,
                     readable = TRUE,)#,minGSSize = 20,maxGSSize = 500)
dotplot(egoHoxc9,showCategory=15)

goplot(egoHoxc9,showCategory=3)




egoHoxc9BP <- enrichGO(gene = testIDs, 
                       keyType = "GID", 
                       OrgDb = org.Aaegypti.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.9,
                       readable = TRUE,minGSSize = 20,maxGSSize = 80)
dotplot(egoHoxc9BP,showCategory=15)

goplot(egoHoxc9BP,showCategory=3)

library(org.Aaegypti.eg.db)

#CIGenes.txt
peaks_CI=read.table("CIGenes.txt")

egoCI <- enrichGO(gene = unlist(peaks_CI), 
                  keyType = "GID", 
                  OrgDb = org.Aaegypti.eg.db, 
                  ont = "MF", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 1,
                  readable = TRUE)#,minGSSize = 10,maxGSSize = 100)
enrichplot::dotplot(egoCI)

CIRes=egoCI@result

nrow(peaks_CI)
sum(unique(peaks_CI$V1))

egoCIBP <- enrichGO(gene = unlist(peaks_CI), 
                    keyType = "GID", 
                    OrgDb = org.Aaegypti.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 1,
                    readable = TRUE)#,minGSSize = 10,maxGSSize = 100)
dotplot(egoCIBP)

CIResBP=egoCIBP@result

edox2 <- pairwise_termsim(egoCIBP)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

goplot(egoCIBP)


#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html




egoHoxcBPDef <- enrichGO(gene = testIDs, 
                         keyType = "GID", 
                         OrgDb = org.Aaegypti.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.9,
                         readable = TRUE)#,minGSSize = 20,maxGSSize = 80)


p1 <- heatplot(egoHoxcBPDef, showCategory=5)
p2 <- heatplot(egoHoxcBPDef, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

edox2 <- pairwise_termsim(egoHoxcBPDef)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')


edo <- pairwise_termsim(egoHoxcBPDef)
p1 <- emapplot(edo)
#p2 <- emapplot(edo, cex_category=1.5)
p3 <- emapplot(edo, layout="kk")
#p4 <- emapplot(edo, cex_category=1.5,layout="kk")  
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])


library(clusterProfiler)
data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
xx <- pairwise_termsim(xx)                     
p1 <- emapplot(xx)
p2 <- emapplot(xx, legend_n=2) 
p3 <- emapplot(xx, pie="count")
p4 <- emapplot(xx, pie="count", cex_category=1.5, layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])



terms <- egoHoxcBPDef$Description[1:5]
p <- pmcplot(terms, 2010:2024)
p2 <- pmcplot(terms, 2010:2024, proportion=FALSE)
plot_grid(p, p2, ncol=2)




egoRSF1BPSelect

terms <- egoRSF1BPSelect$Description[1:5]
p <- pmcplot(terms, 2010:2024)
p2 <- pmcplot(terms, 2010:2024, proportion=FALSE)
plot_grid(p, p2, ncol=2)



edox2 <- pairwise_termsim(egoRSF1BP)#egoRSF1BPSelect)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

ego <- pairwise_termsim(egoRSF1BPSelect)
p1 <- treeplot(ego, hclust_method = "average")
p2 <- treeplot(edox2, hclust_method = "average")

peaks_RFX7=read.table("RFX7Genes.txt")
peaks_RFX7=unique(peaks_RFX7)
peaks_RFX7$
  
  egoRFX7 <- enrichGO(gene = peaks_RFX7$V1, 
                      keyType = "GID", 
                      OrgDb = org.Aaegypti.eg.db, 
                      ont = "CC", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 1,
                      readable = TRUE)#,minGSSize = 1,maxGSSize = 100)
dotplot(egoRFX7)



compareChIPGSEA=function(fileList,orgData,txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",rankMetric = "pval",goCat="All",qval=0.5,plotting=FALSE,sampNames=c(),promoterOnly=TRUE,enhancerOnly=FALSE,dists=c(-2000,2000),keycol="GID"){
  #myRankedGenes=vector("list", length = length(fileList))
  tempFile="tempDiffbind.tsv"
  testtx=makeTxDbFromGFF(txFile)
  count=0
  chipgseahelper=function(diffbindFile){
    
    myDiffbind=read.csv(diffbindFile)
    slimDiffbind=myDiffbind[c("chr","start","end")]
    
    slimDiffbind["name"]=myDiffbind[[1]]
    slimDiffbind["score"]=0
    slimDiffbind["strand"]="."
    slimDiffbind["signalValue"]=myDiffbind$Fold
    slimDiffbind["pValue"]=-log10(myDiffbind$p.value)
    slimDiffbind["qValue"]=-log10(myDiffbind$FDR)
    slimDiffbind["peak"]=-1
    
    write.table(slimDiffbind,tempFile,sep="\t",row.names=FALSE)
    
    
    myPeaks=annotatePeak(tempFile,tssRegion = dists,TxDb=testtx,verbose=FALSE)
    
    ex_annot = as.data.frame(myPeaks)#myPeaks@anno
    print("The number of rows pre-filtering is")
    print(nrow(ex_annot))
    if(promoterOnly){
      ex_annot=ex_annot[abs(ex_annot$distanceToTSS) <= dists[2],]
    }
    #Best way to incorporate distances?
    #Should this override promoter if TRUE since non-default?
    else if(enhancerOnly){
      ex_annot=ex_annot[abs(ex_annot$distanceToTSS) >= 50000,]
      ex_annot=ex_annot[abs(ex_annot$distanceToTSS) <= 200000,]
      ex_annot=ex_annot[ex_annot$annotation=="Distal Intergenic",]
    }
    print("The number of rows after filtering is")
    print(nrow(ex_annot))
    #Initially sorted by -log10(pvalue) to keep most significant gene for category, but could change for other methods 
    ex_annot$pValue
    ex_annot=ex_annot %>%
      arrange(desc(pValue))
    
    unique_annot=ex_annot[!duplicated(ex_annot$geneId),]
    
    print("After adjusting for uniqueness")
    print(nrow(unique_annot))
    
    if(rankMetric=="pval"){
      unique_annot=arrange(unique_annot,desc(unique_annot$pValue*sign(unique_annot$signalValue)))
      
      myranked=unique_annot$pValue*sign(unique_annot$signalValue)
    }
    
    if(rankMetric=="foldchange"){
      unique_annot=arrange(unique_annot,desc(unique_annot$signalValue))
      
      myranked=unique_annot$signalValue
    }
    
    if(rankMetric=="both"){
      unique_annot=arrange(unique_annot,desc(unique_annot$pValue*unique_annot$signalValue))
      myranked=unique_annot$pValue*unique_annot$signalValue
    }
    names(myranked)=unique_annot$geneId
    return(myranked)
    #myRankedGenes[count]=myranked
    #count=count+1
  }
  
  #names(my)
  #Significantly better results with duplicate peaks
  # ego <- gseGO(geneList     = myRankedGenes,
  #              OrgDb        = orgData,
  #              keyType="GID",
  #              ont          = goCat,
  #              minGSSize    = minSize,
  #              maxGSSize    = maxSize,
  #              pvalueCutoff = 0.9,
  #              verbose      = TRUE)
  #myRankedGenes=list(myRankedGenes)
  #names(myRankedGenes)=sampNames
  #print(myRankedGenes)
  pre_genes = lapply(fileList,chipgseahelper)#function(i) as.data.frame(i)$geneId)
  
  if(length(sampNames)==0){
    sampNames=paste0("Sample_", 1:length(fileList))
  }
  
  names(pre_genes)=sampNames
  set.seed(2025)
  compGO <- compareCluster(geneCluster = pre_genes,#genes, 
                           #universe=backgroundids,
                           fun = "gseGO",
                           keyType = "GID", 
                           OrgDb = orgData, 
                           ont = goCat, 
                           pvalueCutoff = qval,
                           pAdjustMethod = "BH")#, readable=TRUE)
  if(plotting){
    dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis") 
  }
  if(file.exists(tempFile)) {
    file.remove(tempFile)
  }
  return(compGO)
}



RVFVd1_H3K27ac_vs_BF-DiffBindout.csv

library(org.Aaegypti.eg.db)






peaks_Irf6=read.table("Irf6Genes.txt")

egoIrf6 <- enrichGO(gene = unlist(peaks_Irf6), 
                    keyType = "GID", 
                    OrgDb = org.Aaegypti.eg.db, 
                    ont = "MF", 
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05,
                    readable = TRUE)#,minGSSize = 10,maxGSSize = 100)
enrichplot::dotplot(egoIrf6)

Irf6Res=egoIrf6@result

edox2 <- pairwise_termsim(egoIrf6)
p1 <- treeplot(edox2)
p1

#nrow(peaks_CI)
sum(unique(peaks_Irf6$V1))

egoIrf6BP <- enrichGO(gene = unlist(peaks_Irf6), 
                      keyType = "GID", 
                      OrgDb = org.Aaegypti.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 1,
                      readable = TRUE)#,minGSSize = 10,maxGSSize = 100)
dotplot(egoIrf6BP)

Irf6ResBP=egoIrf6BP@result

egoIrf6CC <- enrichGO(gene = unlist(peaks_Irf6), 
                      keyType = "GID", 
                      OrgDb = org.Aaegypti.eg.db, 
                      ont = "CC", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 1,
                      readable = TRUE)#,minGSSize = 10,maxGSSize = 100)
dotplot(egoIrf6CC)

Irf6ResCC=egoIrf6CC@result

edox2 <- pairwise_termsim(egoIrf6BP)
p1 <- treeplot(edox2)
p1
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

goplot(egoCIBP)












ChIPGSEA=function(peakFile,orgData,rankMetric = "pval",txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",goCat="All",qval=0.5,promoterOnly=TRUE,enhancerOnly=FALSE,minSize=15,maxSize=500,plotting=TRUE,dists=c(-2000,2000)){
  
  tempFile="tempDiffbind.tsv"
  
  #testtx=makeTxDbFromGFF(txFile)
  
  myDiffbind=read.csv(peakFile)
  slimDiffbind=myDiffbind[c("chr","start","end")]
  
  slimDiffbind["name"]=myDiffbind[[1]]
  slimDiffbind["score"]=0
  slimDiffbind["peakStrand"]="."
  slimDiffbind["signalValue"]=myDiffbind$Fold
  slimDiffbind["pValue"]=-log10(myDiffbind$p.value)
  slimDiffbind["qValue"]=-log10(myDiffbind$FDR)
  slimDiffbind["peak"]=-1
  
  write.table(slimDiffbind,tempFile,sep="\t",row.names=FALSE)
  
  testtx=makeTxDbFromGFF(txFile)
  myPeaks=annotatePeak(tempFile,tssRegion = dists,TxDb=testtx,verbose=FALSE)
  
  print("PastPeakAnnotation!")
  
  ex_annot = as.data.frame(myPeaks)#myPeaks@anno
  print("The number of rows pre-filtering is")
  print(nrow(ex_annot))
  if(promoterOnly){
    ex_annot=ex_annot[abs(ex_annot$distanceToTSS) <= dists[2],]
  }
  #Best way to incorporate distances?
  #Should this override promoter if TRUE since non-default?
  else if(enhancerOnly){
    ex_annot=ex_annot[abs(ex_annot$distanceToTSS) >= 50000,]
    ex_annot=ex_annot[abs(ex_annot$distanceToTSS) <= 200000,]
    ex_annot=ex_annot[ex_annot$annotation=="Distal Intergenic",]
  }
  print("The number of rows after filtering is")
  print(nrow(ex_annot))
  #Initially sorted by -log10(pvalue) to keep most significant gene for category, but could change for other methods 
  #ex_annot$pValue
  ex_annot=ex_annot %>%
    arrange(desc(pValue))
  #return(ex_annot)
  unique_annot=ex_annot[!duplicated(ex_annot$geneId),]
  
  print("After adjusting for uniqueness")
  print(nrow(unique_annot))
  
  if(rankMetric=="pval"){
    unique_annot=arrange(unique_annot,desc(unique_annot$pValue*sign(unique_annot$signalValue)))
    
    myranked=unique_annot$pValue*sign(unique_annot$signalValue)
  }
  
  if(rankMetric=="foldchange"){
    unique_annot=arrange(unique_annot,desc(unique_annot$signalValue))
    
    myranked=unique_annot$signalValue
  }
  
  if(rankMetric=="both"){
    unique_annot=arrange(unique_annot,desc(unique_annot$pValue*unique_annot$signalValue))
    myranked=unique_annot$pValue*unique_annot$signalValue
  }
  print("Ready to GO!")
  
  names(myranked)=unique_annot$geneId
  
  #Significantly better results with duplicate peaks
  set.seed(2025)
  ego <- gseGO(geneList     = myranked,
               OrgDb        = orgData,
               keyType="GID",
               ont          = goCat,
               minGSSize    = minSize,
               maxGSSize    = maxSize,
               pvalueCutoff = 0.9,
               verbose      = TRUE)
  if(plotting){
    dotplot(ego) 
  }
  
  if(file.exists(tempFile)) {
    file.remove(tempFile)
  }
  return(ego)
}

#Example Usage
#The input file is the output of diffbind, written to a csv

#Using only promoter regions
mychip=ChIPGSEA("BF_H3K27ac_d7_vs_SF_H3K27ac-DiffBind-out-NEW.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",rankMetric = "both")

mychipMF=ChIPGSEA("BF_H3K27ac_d7_vs_SF_H3K27ac-DiffBind-out-NEW.csv",orgData = org.Aaegypti.eg.db,goCat = "MF",rankMetric = "pval")
dotplot(mychipMF)

dotplot(mychipMF,split = ".sign") + facet_grid(~.sign)

goplot(mychipMF,split = ".sign") + facet_grid(~.sign)









dotplot(mychip)
bpres=mychip@result

#unique_annot=arrange(mychip,desc(mychip$signalValue))
#duplicated(colnames(mychip))
#Had strand name in multiple
mychip$strand

#Using only putative enhancer regions
bpchip=ChIPGSEA("RVFVd1_H3K27ac_vs_BF-DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE)#,enhancerOnly = TRUE)
#BF_H3K27ac_d7_vs_SF_H3K27ac-DiffBind-out-NEW.csv
dotplot(bpchip)


d1diffbind=read.csv("RVFV_H3K9me_d1_vs_BF_p10-DiffBindout.csv")

d3diffbind=read.csv("RVFV_H3K9me_d3_vs_BF-DiffBindout.csv")


d7diffbind=read.csv("RVFV_H3K9me_vs_BF-keepdup-FINAL-DiffBind-out.csv")


d1MeBP=ChIPGSEA("RVFV_H3K9me_d1_vs_BF_p10-DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = TRUE)#,enhancerOnly = TRUE)
dotplot(d1MeBP)

d3MeBP=ChIPGSEA("RVFV_H3K9me_d3_vs_BF-DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = TRUE)#,enhancerOnly = TRUE)
dotplot(d3MeBP)
edox2 <- pairwise_termsim(d3MeBP)
p1 <- treeplot(edox2)
p1


d7MeBP=ChIPGSEA("RVFV_H3K9me_vs_BF-keepdup-FINAL-DiffBind-out.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = TRUE)#,enhancerOnly = TRUE)
dotplot(d7MeBP)

d7MeBPRes=d7MeBP@result

edox2 <- pairwise_termsim(d7MeBP)
p1 <- treeplot(edox2)
p1

d7MeMF=ChIPGSEA("RVFV_H3K9me_vs_BF-keepdup-FINAL-DiffBind-out.csv",orgData = org.Aaegypti.eg.db,goCat = "MF",promoterOnly = TRUE)#,enhancerOnly = TRUE)
dotplot(d7MeMF)


edox2 <- pairwise_termsim(d7MeMF)
p1 <- treeplot(edox2)
p1


d7MeBPSel=ChIPGSEA("RVFV_H3K9me_vs_BF-keepdup-FINAL-DiffBind-out.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = TRUE,maxSize = 80)#,enhancerOnly = TRUE)
dotplot(d7MeBPSel)


edox2 <- pairwise_termsim(d7MeBPSel)
p1 <- treeplot(edox2)
p1



p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')





peaks_ISL2=read.table("ISL2Genes.txt")

egoISL2 <- enrichGO(gene = unlist(peaks_ISL2), 
                    keyType = "GID", 
                    OrgDb = org.Aaegypti.eg.db, 
                    ont = "MF", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 1,
                    readable = TRUE)#,minGSSize = 10,maxGSSize = 100)
dotplot(egoISL2)


egoISL2CC <- enrichGO(gene = unlist(peaks_ISL2), 
                      keyType = "GID", 
                      OrgDb = org.Aaegypti.eg.db, 
                      ont = "CC", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 1,
                      readable = TRUE)#,minGSSize = 10,maxGSSize = 100)
dotplot(egoISL2CC)

"ISL2Genes.txt"
edox2 <- pairwise_termsim(egoISL2CC)
p1 <- treeplot(edox2)
p1
#Putative methylation



d1MeBP=ChIPGSEA("RVFV_H3K9me_d1_vs_BF_p10-DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)
dotplot(d1MeBP)

d3MeBP=ChIPGSEA("RVFV_H3K9me_d3_vs_BF-DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "All",promoterOnly = FALSE,enhancerOnly = TRUE)
dotplot(d3MeBP)
edox2 <- pairwise_termsim(d3MeBP)
p1 <- treeplot(edox2)
p1


d7MeBP=ChIPGSEA("RVFV_H3K9me_vs_BF-keepdup-FINAL-DiffBind-out.csv",orgData = org.Aaegypti.eg.db,goCat = "All",promoterOnly = FALSE,enhancerOnly = TRUE)
#dotplot(d7MeBP)

d7MeBPRes=d7MeBP@result

edox2 <- pairwise_termsim(d7MeBP)
p1 <- treeplot(edox2)
p1

d7MeMF=ChIPGSEA("RVFV_H3K9me_vs_BF-keepdup-FINAL-DiffBind-out.csv",orgData = org.Aaegypti.eg.db,goCat = "MF",promoterOnly = TRUE)#,enhancerOnly = TRUE)
dotplot(d7MeMF)


edox2 <- pairwise_termsim(d7MeMF)
p1 <- treeplot(edox2)
p1


d7MeBPSel=ChIPGSEA("RVFV_H3K9me_vs_BF-keepdup-FINAL-DiffBind-out.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = TRUE,maxSize = 80)#,enhancerOnly = TRUE)
dotplot(d7MeBPSel)


edox2 <- pairwise_termsim(d7MeBPSel)
p1 <- treeplot(edox2)
p1


getGenesRegion=function(peakFile,orgData,rankMetric = "pval",txFile="VectorBase-68_AaegyptiLVP_AGWG.gff",goCat="All",qval=0.5,promoterOnly=TRUE,enhancerOnly=FALSE,minSize=15,maxSize=500,plotting=TRUE,dists=c(-2000,2000)){
  
  tempFile="tempDiffbind.tsv"
  
  #testtx=makeTxDbFromGFF(txFile)
  
  myDiffbind=read.csv(peakFile)
  slimDiffbind=myDiffbind[c("chr","start","end")]
  
  slimDiffbind["name"]=myDiffbind[[1]]
  slimDiffbind["score"]=0
  slimDiffbind["peakStrand"]="."
  slimDiffbind["signalValue"]=myDiffbind$Fold
  slimDiffbind["pValue"]=myDiffbind$p.value#-log10(myDiffbind$p.value)
  slimDiffbind["qValue"]=myDiffbind$FDR#-log10(myDiffbind$FDR)
  slimDiffbind["peak"]=-1
  
  write.table(slimDiffbind,tempFile,sep="\t",row.names=FALSE)
  
  testtx=makeTxDbFromGFF(txFile)
  myPeaks=annotatePeak(tempFile,tssRegion = dists,TxDb=testtx,verbose=FALSE)
  
  print("PastPeakAnnotation!")
  
  ex_annot = as.data.frame(myPeaks)#myPeaks@anno
  print("The number of rows pre-filtering is")
  print(nrow(ex_annot))
  if(promoterOnly){
    ex_annot=ex_annot[abs(ex_annot$distanceToTSS) <= dists[2],]
  }
  #Best way to incorporate distances?
  #Should this override promoter if TRUE since non-default?
  else if(enhancerOnly){
    ex_annot=ex_annot[abs(ex_annot$distanceToTSS) >= 50000,]
    ex_annot=ex_annot[abs(ex_annot$distanceToTSS) <= 200000,]
    ex_annot=ex_annot[ex_annot$annotation=="Distal Intergenic",]
  }
  print("The number of rows after filtering is")
  print(nrow(ex_annot))
  #Initially sorted by -log10(pvalue) to keep most significant gene for category, but could change for other methods 
  #ex_annot$pValue
  ex_annot=ex_annot %>%
    arrange(desc(pValue))
  #return(ex_annot)
  unique_annot=ex_annot[!duplicated(ex_annot$geneId),]
  
  print("After adjusting for uniqueness")
  print(nrow(unique_annot))
  
  if(rankMetric=="pval"){
    unique_annot=arrange(unique_annot,desc(unique_annot$pValue*sign(unique_annot$signalValue)))
    
    myranked=unique_annot$pValue*sign(unique_annot$signalValue)
  }
  
  if(rankMetric=="foldchange"){
    unique_annot=arrange(unique_annot,desc(unique_annot$signalValue))
    
    myranked=unique_annot$signalValue
  }
  
  if(rankMetric=="both"){
    unique_annot=arrange(unique_annot,desc(unique_annot$pValue*unique_annot$signalValue))
    myranked=unique_annot$pValue*unique_annot$signalValue
  }
  print("Ready to GO!")
  return(unique_annot)
  names(myranked)=unique_annot$geneId
  
  #Significantly better results with duplicate peaks
  set.seed(2025)
  ego <- gseGO(geneList     = myranked,
               OrgDb        = orgData,
               keyType="GID",
               ont          = goCat,
               minGSSize    = minSize,
               maxGSSize    = maxSize,
               pvalueCutoff = 0.9,
               verbose      = TRUE)
  if(plotting){
    dotplot(ego) 
  }
  
  if(file.exists(tempFile)) {
    file.remove(tempFile)
  }
  return(ego)
}

d7meEnhancers=getGenesRegion("RVFV_H3K9me_vs_BF-keepdup-FINAL-DiffBind-out.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)

d1meEnhancers=getGenesRegion("RVFV_H3K9me_d1_vs_BF_p10-DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)

d3meEnhancers=getGenesRegion("RVFV_H3K9me_d3_vs_BF-DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)


write.csv(d7meEnhancers,"RVFVvBF_D7_Me_DiffbindEnhancer.csv")
write.csv(d3meEnhancers,"RVFVvBF_D3_Me_DiffbindEnhancer.csv")
write.csv(d1meEnhancers,"RVFVvBF_D1_Me_DiffbindEnhancer.csv")




RVFVd1_H3K27ac_vs_BF_DiffBindout.csv
RVFVd3_H3K27ac_vs_BF_DiffBindout.csv
RVFV_H3K27AcvsBF_H3K27Ac.csv


d7acEnhancers=getGenesRegion("RVFV_H3K27AcvsBF_H3K27Ac.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)

d1acEnhancers=getGenesRegion("RVFVd1_H3K27ac_vs_BF_DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)

d3acEnhancers=getGenesRegion("RVFVd3_H3K27ac_vs_BF_DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)




write.csv(d7acEnhancers,"RVFVvBF_D7_Ac_DiffbindEnhancer.csv")
write.csv(d3acEnhancers,"RVFVvBF_D3_Ac_DiffbindEnhancer.csv")
write.csv(d1acEnhancers,"RVFVvBF_D1_Ac_DiffbindEnhancer.csv")






write.csv(d7meEnhancers[which(d7meEnhancers$qValue<0.1),],"RVFVvBF_D7_Me_FDR_0.1_DiffbindEnhancer.csv")
write.csv(d3meEnhancers[which(d3meEnhancers$qValue<0.1),],"RVFVvBF_D3_Me_FDR_0.1_DiffbindEnhancer.csv")
write.csv(d1meEnhancers[which(d1meEnhancers$qValue<0.1),],"RVFVvBF_D1_Me_FDR_0.1_DiffbindEnhancer.csv")

write.csv(d7acEnhancers[which(d7acEnhancers$qValue<0.1),],"RVFVvBF_D7_FDR_0.1_AcDiffbindEnhancer.csv")
write.csv(d3acEnhancers[which(d3acEnhancers$qValue<0.1),],"RVFVvBF_D3_FDR_0.1_DiffbindEnhancer.csv")
write.csv(d1acEnhancers[which(d1acEnhancers$qValue<0.1),],"RVFVvBF_D1_FDR_0.1_DiffbindEnhancer.csv")




#BF_H3K9Me3vsSF_H3K9Me3_Day7.csv
#BF_H3K9Me3vsSF_H3K9Me3_Day3.csv

d7meEnhancers=getGenesRegion("BF_H3K9Me_d7_vs_SF_H3K9Me-DiffBind-out.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)

d3meEnhancers=getGenesRegion("BF_H3K9Me3vsSF_H3K9Me3_Day3.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)

d1meEnhancers=getGenesRegion("BF_H3K9me_d1_vs_SF-keepdup-FINAL-DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)




#RVFVd1_H3K27ac_vs_BF_DiffBindout.csv
#RVFVd3_H3K27ac_vs_BF_DiffBindout.csv
#RVFV_H3K27AcvsBF_H3K27Ac.csv
#BF_H3K27ac_d7_vs_SF_H3K27ac-DiffBind-out-NEW.csv

d7acEnhancers=getGenesRegion("BF_H3K27ac_d7_vs_SF_H3K27ac-DiffBind-out.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)

d3acEnhancers=getGenesRegion("BF_d3_H3K27ac_vs_SF_H3K27ac-DIffBind-out.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)

d1acEnhancers=getGenesRegion("BFd1_H3K27ac_vs_SF-kd-FINAL-DiffBindout.csv",orgData = org.Aaegypti.eg.db,goCat = "BP",promoterOnly = FALSE,enhancerOnly = TRUE)

#problem=read.csv("BFd1_H3K27ac_vs_SF-kd-FINAL-DiffBindout.csv")
#ok=read.csv("BF_d3_H3K27ac_vs_SF_H3K27ac-DIffBind-out.csv")


write.csv(d7meEnhancers,"BFvSF_D7_Me_DiffbindEnhancer.csv")
write.csv(d3meEnhancers,"BFvSF_D3_Me_DiffbindEnhancer.csv")
write.csv(d1meEnhancers,"BFvSF_D1_Me_DiffbindEnhancer.csv")



write.csv(d7acEnhancers,"BFvSF_D7_Ac_DiffbindEnhancer.csv")
write.csv(d3acEnhancers,"BFvSF_D3_Ac_DiffbindEnhancer.csv")
write.csv(d1acEnhancers,"BFvSF_D1_Ac_DiffbindEnhancer.csv")


write.csv(d7meEnhancers[which(d7meEnhancers$qValue<0.1),],"BFvSF_D7_Me_FDR_0.1_DiffbindEnhancer.csv")
write.csv(d3meEnhancers[which(d3meEnhancers$qValue<0.1),],"BFvSF_D3_Me_FDR_0.1_DiffbindEnhancer.csv")
write.csv(d1meEnhancers[which(d1meEnhancers$qValue<0.1),],"BFvSF_D1_Me_FDR_0.1_DiffbindEnhancer.csv")

write.csv(d7acEnhancers[which(d7acEnhancers$qValue<0.1),],"BFvSF_D7_FDR_0.1_AcDiffbindEnhancer.csv")
write.csv(d3acEnhancers[which(d3acEnhancers$qValue<0.1),],"BFvSF_D3_FDR_0.1_DiffbindEnhancer.csv")
write.csv(d1acEnhancers[which(d1acEnhancers$qValue<0.1),],"BFvSF_D1_FDR_0.1_DiffbindEnhancer.csv")






RVFVd7acEnhancers=read.csv("RVFVvBF_D7_Ac_DiffbindEnhancer.csv")
RVFVd3acEnhancers=read.csv("RVFVvBF_D3_Ac_DiffbindEnhancer.csv")
RVFVd1acEnhancers=read.csv("RVFVvBF_D1_Ac_DiffbindEnhancer.csv")



RVFVd7meEnhancers=read.csv("BFvSF_D7_Me_DiffbindEnhancer.csv")
RVFVd3meEnhancers=read.csv("BFvSF_D3_Me_DiffbindEnhancer.csv")
RVFVd1meEnhancers=read.csv("BFvSF_D1_Me_DiffbindEnhancer.csv")



SFd7acEnhancers=read.csv("BFvSF_D7_Ac_DiffbindEnhancer.csv")
SFd3acEnhancers=read.csv("BFvSF_D3_Ac_DiffbindEnhancer.csv")
SFd1acEnhancers=read.csv("BFvSF_D1_Ac_DiffbindEnhancer.csv")


SFd7meEnhancers=read.csv("RVFVvBF_D7_Me_DiffbindEnhancer.csv")
SFd3meEnhancers=read.csv("RVFVvBF_D3_Me_DiffbindEnhancer.csv")
SFd1meEnhancers=read.csv("RVFVvBF_D1_Me_DiffbindEnhancer.csv")

nrow(SFd7meEnhancers[which(SFd7meEnhancers$qValue<0.1),])

RVFVd7RNA=read.csv("RVFVvBF_3d_DESeqresults_min_10_all-annot.csv")
RVFVd3RNA=read.csv("RVFVvBF_DAY1_DESeqresults_min_10_all-annot.csv")
RVFVd1RNA=read.csv("RVFVvBF7d_DESEQ2p10.csv")

SFd7RNA=read.csv("DAY7_BFvSFdresults.csv")
SFd3RNA=read.csv("DAY3_BFvSFdresultsp10_ms-PEannotated.csv")
SFd1RNA=read.csv("DAY1_BFvSFdresultsp10_ms-PEannotated.csv")

acEn=list(RVFVd1acEnhancers,RVFVd3acEnhancers,RVFVd7acEnhancers,SFd1acEnhancers,SFd3acEnhancers,SFd7acEnhancers)

meEn=list(RVFVd1meEnhancers,RVFVd3meEnhancers,RVFVd7meEnhancers,SFd1meEnhancers,SFd3meEnhancers,SFd7meEnhancers)

rnas=list(RVFVd1RNA,RVFVd3RNA,RVFVd7RNA,SFd1RNA,SFd3RNA,SFd7RNA)

df <- data.frame(matrix(ncol = 6, nrow = 12))
dfcols <- c("Total Differentially Bound Peaks In Putative Enhancer Regions", "Prop Peaks enriched", "Prop Peaks depleted","# Peaks depleted in RVFV exposed","# DEGs assoc w/ peaks","total DEGs")
colnames(df) <- dfcols

#RVFVd3acEnhancers$qValue
SFd1acEnhancers$geneId%in%SFd1RNA$X
#SFd1acEnhancers
#SFd1RNA
#acEn$

#tempAc=acEn[[1]][which(acEn[[1]]$qValue <0.1),]
for(i in seq(6)){
  tempAc=acEn[[i]][which(acEn[[i]]$qValue <0.1),]
  tempMe=meEn[[i]][which(meEn[[i]]$qValue <0.1),]
  tempRNA=rnas[[i]][which(rnas[[i]]$padj <0.1),]
  #Total
  print(nrow(acEn[[i]]))
  print(nrow(meEn[[i]]))
  df[i*2-1,1]=nrow(tempAc)
  df[i*2,1]=nrow(tempMe)
  #Percent Enriched
  df[i*2-1,2]=nrow(tempAc[which(tempAc$signalValue>0),])/nrow(tempAc)
  df[i*2,2]=nrow(tempMe[which(tempMe$signalValue>0),])/nrow(tempMe)
  #Percent Depleted
  df[i*2-1,3]=nrow(tempAc[which(tempAc$signalValue<0),])/nrow(tempAc)
  df[i*2,3]=nrow(tempMe[which(tempMe$signalValue<0),])/nrow(tempMe)
  #NumberDepleted
  df[i*2-1,4]=nrow(tempAc[which(tempAc$signalValue<0),])
  df[i*2,4]=nrow(tempMe[which(tempMe$signalValue<0),])
  
  #ToDEGs
  df[i*2-1,5]=sum(tempAc$geneId%in%tempRNA$geneID)
  df[i*2,5]=sum(tempMe$geneId%in%tempRNA$geneID)
  
  df[i*2-1,6]=nrow(tempRNA)
  df[i*2,6]=nrow(tempRNA)
  #=nrow()
  #tempMe=
  
}

is.num <- sapply(df, is.numeric)
df[is.num] <- lapply(df[is.num], round, 3)

write.csv(df,"EnhancerMetadataCalc.csv")
