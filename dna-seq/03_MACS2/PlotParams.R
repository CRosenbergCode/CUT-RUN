library(readr)
library(ggplot2)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Read in tsv created by bash script
bf_ago2_d7_df=read_tsv("BF_Ago2_Merged_D7ParamSummary.tsv")#"BF_Ac_1_ParamSummaryLengthAdjusted.tsv")

sf_ago2_d7_df=read_tsv("SF_Ago2_Merged_D7_ParamSummary.tsv")

#Get median number of peaks
median(bf_ago2_d7_df$TotalPeaks)
median(sf_ago2_d7_df$TotalPeaks)

#Plot the relationship between two columns
ggplot(data=bf_ago2_d7_df)+geom_point(aes(x=Slocal, y=AedesGenesProm1500.gtf))

#Get the median value of one column based on a given value for another
median(bf_ago2_d7_df[bf_ago2_d7_df$QValue==0.05,]$`Median Peak Width`)

#Set variables to be treated as factors, rather than numbers
#VERY IMPORTANT FOR THE REGRESSION, DO NOT SKIP THIS
bf_ago2_d7_df$QValue=factor(bf_ago2_d7_df$QValue)
bf_ago2_d7_df$D=factor(bf_ago2_d7_df$D)
bf_ago2_d7_df$Llocal=factor(bf_ago2_d7_df$Llocal)
bf_ago2_d7_df$Slocal=factor(bf_ago2_d7_df$Slocal)

#Create a linear model with the number of peaks as the outcome and all other variables as predictors
BF_Ago2_Mod.peaks <- lm(TotalPeaks ~ QValue + D + Slocal + Llocal, data = bf_ago2_d7_df)

#Create a linear model with a single predictor
BF_Ago2_Mod.10k <- lm(AedesGenesProm10000.gtf ~ QValue, data = bf_ago2_d7_df)

#Get summaries for a model, including the coefficients and variance
summary(BF_Ago2_Mod.10k)






ggplot(data=sf_ago2_d7_df)+geom_point(aes(x=Slocal, y=AedesGenesProm1500.gtf))

median(sf_ago2_d7_df[sf_ago2_d7_df$QValue==0.05,]$`Median Peak Width`)
median(sf_ago2_d7_df[sf_ago2_d7_df$QValue==0.00001,]$`Median Peak Width`)
median(sf_ago2_d7_df[sf_ago2_d7_df$QValue==0.05,]$TotalPeaks)
median(sf_ago2_d7_df[sf_ago2_d7_df$QValue==0.005,]$`Total Peaks`)
median(sf_ago2_d7_df[sf_ago2_d7_df$QValue==0.00001,]$`Total Peaks`)

sf_ago2_d7_df$QValue=factor(sf_ago2_d7_df$QValue)
sf_ago2_d7_df$D=factor(sf_ago2_d7_df$D)
sf_ago2_d7_df$Llocal=factor(sf_ago2_d7_df$Llocal)
sf_ago2_d7_df$Slocal=factor(sf_ago2_d7_df$Slocal)

sf_Ago2_Mod.peaks <- lm(TotalPeaks ~ QValue + D + Slocal + Llocal, data = sf_ago2_d7_df)

sf_Ago2_Mod.10k <- lm(AedesGenesProm10000.gtf ~ QValue, data = sf_ago2_d7_df)

sf_Ago2_Mod.10k <- lm(AedesGenesProm1000.gtf ~ QValue, data = sf_ago2_d7_df)



sf_Ago2_Mod.10k <- lm(MedianPeakWidth ~ QValue, data = sf_ago2_d7_df)

sf_Ago2_Mod.10k <- lm(AedesGenesProm10000.gtf ~ D, data = sf_ago2_d7_df)
summary(sf_Ago2_Mod.10k)

sf_Ago2_Mod.10k <- lm(MedianPeakWidth ~ D, data = sf_ago2_d7_df)
summary(sf_Ago2_Mod.10k)


summary(sf_Ago2_Mod.10k)

summary(sf_Ago2_Mod.peaks)


mydata$rank <- factor(mydata$rank)
mylogit <- glm(admit ~ gre + gpa + rank, data = mydata, family = "binomial")