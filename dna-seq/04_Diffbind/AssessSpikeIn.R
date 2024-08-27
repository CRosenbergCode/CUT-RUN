yeastNorm = c(2520,7053,2290,10253,1140,33174,2914,7089,6326,884,9408,3564,61606)
coliNorm = c(713,617,411,517,223,2739,322,278,633,77,372,122,4910)
yeastNormRep2 = c(7053,10253,33174,6326,9408,61606)
coliNormRep2 = c(617,517,2739,633,372,4910)

length(yeastNorm)
length(coliNorm)
cor.test(yeastNorm,coliNorm,method=c("pearson"))#, "kendall", "spearman"))

library(ggplot2)

library("ggpubr")
# mpg
ggqqplot(yeastNorm, ylab = "S. cerevisiae")
# wt
ggqqplot(coliNorm, ylab = "E. coli")

logColi = log(coliNorm)
logYeast = log(yeastNorm)

test = data.frame(yeastNorm,coliNorm,logYeast,logColi)

logColi2 = log(coliNormRep2)
logYeast2 = log(yeastNormRep2)

test2 = data.frame(yeastNormRep2,coliNormRep2,logYeast2,logColi2)


ggscatter(test, x = "yeastNorm", y = "coliNorm", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Yeast Counts", ylab = "E Coli Counts")# + xlim(scale=log10)

ggscatter(test, x = "logYeast", y = "logColi", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Yeast Counts", ylab = "E Coli Counts") 

ggscatter(test2, x = "yeastNormRep2", y = "coliNormRep2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Yeast Counts", ylab = "E Coli Counts")

ggscatter(test2, x = "logYeast2", y = "logColi2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Yeast Counts", ylab = "E Coli Counts")