library(ggplot2)
library(RColorBrewer)
library(reshape2)

setwd("C:\\Users\\hunte\\Desktop\\ChipData")

sig_genes_annotations <- read.csv("sig_genes_annotations_padj05.csv")
custom_labels <- c("SIG" = "Cell-cell interactions/signaling",
                   "Detox/ReDox" = "ReDox - detoxification",
                   "DIV" = "Diverse function",
                   "MET" = "Metabolism or lipid-associated",
                   "PROT" = "Proteolysis or proteosomal activity",
                   "RRTT" = "Replication, Transcription, Translation, DNA repair",
                   "TRP" = "Transport",
                   "UNK" = "Unknown function",
                   "UNK" = "Unknown function", "UNK" = "Unknown function")
# Modify the ggplot2 code to use scale_fill_manual with custom labels
Subset <- sig_genes_annotations$Functional
summary_data <- aggregate(log2FoldChange ~ Gene.ID + Functional, data = sig_genes_annotations, FUN = mean)
temp_data = summary_data[order(summary_data$Functional,summary_data$log2FoldChange),]
#temp_data$Gene.ID <- factor(temp_data$Gene.ID, levels = temp_data$log2FoldChange)
#melt_data = melt(summary_data)
temp_data$Number = seq(length(temp_data$Gene.ID))
g = ggplot(temp_data, aes(x = Gene.ID, y = log2FoldChange, fill = Functional)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  xlab("Gene ID") +
  ylab("Mean log2 Fold Change") +
  ggtitle("Mean log2 Fold Change by Gene ID and Function") +
  facet_wrap(vars(temp_data$Functional)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = brewer.pal(length(custom_labels), "Set3"),
                    name = "Function",
                    labels = custom_labels) 
  #scale_x_continuous(breaks=temp_data$Number, labels=temp_data$Gene.ID)


g = ggplot(temp_data, aes(x = Number, y = log2FoldChange, fill = Functional)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  xlab("Gene ID") +
  ylab("Mean log2 Fold Change") +
  ggtitle("Mean log2 Fold Change by Gene ID and Function") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = brewer.pal(length(custom_labels), "Set3"),
                  name = "Function",
                  labels = custom_labels) + 
  #annotate("label", x = 1, y = 8, label = "CYT/STR") +
  annotate("label", x = 4, y = 8, label = "Redox") +
  annotate("label", x = 15, y = 8, label = "Diverse") +
  annotate("label", x = 35, y = 8, label = "Metabolism") +
  annotate("label", x = 57, y = 8, label = "Transport") +
  annotate("label", x = 68, y = 8, label = "Unknown") + scale_x_continuous(breaks=temp_data$Number, labels=temp_data$Gene.ID)


ggplot(melt_data, aes(x = Gene.ID, y = log2FoldChange, fill = Functional)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  xlab("Gene ID") +
  ylab("Mean log2 Fold Change") +
  ggtitle("Mean log2 Fold Change by Gene ID and Function") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = brewer.pal(length(custom_labels), "Set3"),
                    name = "Function",
                    labels = custom_labels)