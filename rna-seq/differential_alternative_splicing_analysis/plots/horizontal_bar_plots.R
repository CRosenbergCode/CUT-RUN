######## horizontal bar  ##############

# Define the data
A3SS_data <- data.frame(
  category = c("DIV", "RRTT"),
  proportion = c(4/5, 1/5)
)

A5SS_data <- data.frame(
  category = c("DIV", "RRTT"),
  proportion = c(1/2, 1/2)
)

RI_data <- data.frame(
  category = "PROT",
  proportion = 1
)

SE_data <- data.frame(
  functional_category = c("Diverse function", "Unknown function", "Proteolysis or proteosomal activity", "Replication & (DNA) repair", "Transport"),
  proportion = c(1/6, 2/6, 1/6, 1/6, 1/6)
)

# define custom colors for the categories
SE_colors <- c("Diverse function" = "#A6CEE3",
               "Unknown function" = "#F4CAE4",
               "Proteolysis or proteosomal activity" = "#FDBF6F",
               "Replication & (DNA) repair" = "#CAB2D6",
               "Transport" = "#B3E2CD")

# Create the horizontal bar plots
p1 <- ggplot(A3SS_data,aes(x = "A3SS", y = proportion, fill = category)) + 
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = A3SS_colors) +
  xlab("") +
  ylab("Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 1, face = "bold"),
        legend.position = "none", 
        plot.margin = margin(t = -50, b = -50))

p2 <- ggplot(A5SS_data,aes(x = "A5SS", y = proportion, fill = category)) + 
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = A5SS_colors) +
  xlab("") +
  ylab("Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 1, face = "bold"),
        legend.position = "none", 
        plot.margin = margin(t = -50, b = -50))

p3 <- ggplot(RI_data,aes(x = "RI", y = proportion, fill = category)) + 
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = RI_colors) +
  xlab("") +
  ylab("Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 1, face = "bold"),
        legend.position = "none",
        plot.margin = margin(t = -50, b= -50))

p4 <- ggplot(SE_data,aes(x = "SE", y = proportion, fill = functional_category)) + 
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = SE_colors) +
  xlab("") +
  ylab("Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 1, face = "bold"),
        plot.margin = margin(t = -50, b = -50)) + 
  labs(fill = "Gene Functional Categories")

# combine plots
combined_plot <- p1 / p2 / p3 / p4 + plot_layout(guides = "collect", heights = unit(rep(0.8, 4), "inches"))
print(combined_plot)
