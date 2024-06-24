###### colored bar graph #############

custom_labels <- c("CYT/STR" = "CYT/STR - cytoskeletal/structural",
                   "IMM" = "IMM - anti-microbial peptide", 
                   "Detox/ReDox" = "Detox/ReDox - detoxification", 
                   "DIV" = "DIV - diverse function", 
                   "MET" = "MET- metabolism or lipid-associated",
                   "PROT" = "PROT - proteolysis or proteosomal activity",
                   "RRTT" = "RRTT - replication & (DNA) repair",
                   "SIG" = "SIG- signal transduction",
                   "TRP" = "TRP - transport",
                   "UNK" = "UNK - unknown function")

# Define custom colors for each label
custom_colors <- c("CYT/STR" = "#F4A582",
                   "IMM" = "#FFDFD3",
                   "Detox/ReDox" = "#B2DF8A",
                   "DIV" = "#A6CEE3",  
                   "MET" = "#FFD700",  
                   "PROT" = "#D2B48C",
                   "RRTT" = "#FDBF6F", 
                   "SIG" = "#CAB2D6",  
                   "TRP" = "#B3E2CD",  
                   "UNK" = "#F4CAE4") 

# Modify the ggplot2 code to use scale_fill_manual with custom labels
summary_data <- aggregate(log2FoldChange ~ Gene.ID + func, data = sig_genes_annotations, FUN = mean)
ggplot(summary_data, aes(x = Gene.ID, y = log2FoldChange, fill = func)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  xlab("Gene ID") +
  ylab("Mean log2 Fold Change") +
  ggtitle("Mean log2 Fold Change by Gene ID and Function") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = custom_colors, 
                    name = "Function", 
                    labels = custom_labels)
