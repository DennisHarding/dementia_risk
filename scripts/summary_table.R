install.packages("gtsummary")
install.packages("dplyr")
library(dplyr)
library(gtsummary)
library(data.table)

ukb_data <- fread("extended_data_v2.tsv")
df <- tibble(ukb_data) 

#select variables to include in summary table
vars <- c("age_baseline",
          "sex",
          "edu_yrs",
          "gene_apoe",
          "alldem_bin",
          "diabetes",
          "phylo_score"
          )


summary_table <- 
  df %>%
  select(all_of(vars)) %>%
  mutate(
    sex = factor(sex,
                 levels = c(0, 1),
                 labels = c("Female", "Male"))
  ) %>%
  tbl_summary(
    by = sex,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  )

summary_table

