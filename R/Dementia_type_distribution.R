library(ggplot2)

# Data
dementia_data <- data.frame(
  type = c("Alzheimer's",
           "Vascular dementia",
           "Mixed dementia",
           "Lewy bodies",
           "Other"),
  percent = c(60, 15, 10, 10, 5)
)

# Sort first
dementia_data <- dementia_data[order(-dementia_data$percent), ]

# Create ordered labels AFTER sorting
dementia_data$label <- paste0(dementia_data$type, " (", dementia_data$percent, "%)")

# Lock factor order to current row order
dementia_data$label <- factor(dementia_data$label,
                              levels = dementia_data$label)

# Colors in same order
colors <- c("SteelBlue",
            "Firebrick",
            "DarkOliveGreen",
            "DarkOrange",
            "SlateGray")

# Plot
ggplot(dementia_data, aes(x = "", y = percent, fill = label)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colors) +
  theme_void(base_size = 20) +
  labs(title = "Distribution of Dementia Case Types") +
  theme(legend.title = element_blank())
