
#### make separate lists of up and down regulated genes
sig_genes_annotations <- read_excel("/Users/pegaheizad/Desktop/sig_genes_annotations_padj___05.xlsx")
up_genes <- sig_genes_annotations[sig_genes_annotations$log2FoldChange > 0, ]
print(up_genes)

write_xlsx(up_genes, "up_genes.xlsx")

down_genes <- sig_genes_annotations[sig_genes_annotations$log2FoldChange < 0, ]
write_xlsx(down_genes, "down_genes.xlsx")



######PI Charts######################################################

# Define the categories and their counts
up_categories <- c("CYT/STR", "Detox/ReDox", "Detox/ReDox", "DIV", "DIV", "RRTT", 
                   "CYT/STR", "RRTT", "RRTT", "SIG", "TRP", "TRP", "TRP", "TRP", 
                   "TRP", "TRP", "UNK", "RRTT", "RRTT", "UNK")
category_counts_up <- table(up_categories)

# Print the category counts
print(category_counts_up)
# Define the colors for the categories
category_colors <- c("#F4A582", "#B2DF8A", "#A6CEE3", "#FDBF6F", "#CAB2D6", "#B3E2CD", "#F4CAE4")

# Create labels with counts
labels_with_counts <- paste(names(category_counts_up), category_counts_up, sep = " (")

# Add a closing parenthesis to the counts
labels_with_counts <- paste(labels_with_counts, ")", sep = "")

# Plot the pie chart with labels
pie(category_counts_up, 
    main = "Functional gene categories for up regulated genes (padj < 0.05)", 
    col = category_colors, 
    labels = labels_with_counts)




# Define the categories and their counts
down_categories <- c("CYT/STR", "Detox/ReDox", "Detox/ReDox", "Detox/ReDox", "Detox/ReDox", 
                     "IMM", "DIV", "DIV", "IMM", "IMM", "DIV", "DIV", "DIV", "SIG", 
                     "DIV", "DIV", "DIV", "DIV", "IMM", "IMM", "MET", "MET", "MET", 
                     "MET", "MET", "MET", "MET", "MET", "PROT", "PROT", "PROT", 
                     "PROT", "PROT", "PROT", "PROT", "RRTT", "RRTT", "SIG", "SIG", 
                     "SIG", "TRP", "TRP", "TRP", "TRP", "TRP", "TRP", "UNK", "UNK", 
                     "UNK", "UNK")
category_counts_down <- table(down_categories)
# Print the category counts
print(category_counts_down)

# Define the colors for the categories
category_colors <- c("#F4A582", "#B2DF8A", "#A6CEE3", "#FFDFD3", "#FFD700", "#D2B48C", "#FDBF6F", "#CAB2D6", "#B3E2CD", "#F4CAE4")

# Create labels with counts
labels_with_counts <- paste(names(category_counts_down), category_counts_down, sep = " (")

# Add a closing parenthesis to the counts
labels_with_counts <- paste(labels_with_counts, ")", sep = "")

# Plot the pie chart with labels
pie(category_counts_down, 
    main = "Functional gene categories for down regulated genes (padj < 0.05)", 
    col = category_colors, 
    labels = labels_with_counts)

