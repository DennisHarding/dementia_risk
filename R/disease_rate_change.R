library(ggplot2)
library(dplyr)

# Create dataset
df <- data.frame(
  Disease = c("Heart Disease",
              "Stroke",
              "Prostate Cancer",
              "COPD",
              "Dementia"),
  Death_2000 = c(303, 99, 62, 45, 33),
  Death_2021 = c(140, 48, 44, 30, 63)
)

# Calculate percent change
df <- df %>%
  mutate(
    Percent_Change = ((Death_2021 - Death_2000) / Death_2000) * 100
  )

# Create barplot
ggplot(df, aes(x = reorder(Disease, Percent_Change),
               y = Percent_Change,
               fill = Percent_Change > 0)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.9) +
  geom_text(aes(label = paste0(round(Percent_Change, 1), "%")),
            vjust = ifelse(df$Percent_Change > 0, -0.5, 1.5),
            size = 5) +
  scale_fill_manual(values = c("TRUE" = "firebrick",
                               "FALSE" = "steelblue")) +
  labs(
    title = "Percent Change in Deaths per 100,000",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold")
  )
