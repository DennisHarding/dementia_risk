install.packages("survival")
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/ggrepel/ggrepel_0.9.6.tar.gz",
  repos = NULL,
  type = "source"
)
install.packages("survminer")

library(survival)
library(survminer)
library(tibble)
library(dplyr)
library(purrr)

source("scripts/formatting.R")
source("R/plot_spline.R")

types_list <- list(
  AD = "AD",
  VD = "VD",
  NonAD = "NonAD"
)

data_list <- list(
  AD = df_ad,
  VD = df_vd,
  NonAD = df_nonad
)

analysis_tbl <- tibble(
  type = types_list,
  data = data_list
)

# Fit cox models
analysis_tbl <- analysis_tbl %>% 
  mutate(
    fit = map(data, ~ coxph(
      Surv(futime, fail_bin) ~ 
        age_baseline +
        sex + 
        edu +
        gene_apoe +
        smok_ever 
      ,
      data = .x
    ))
  )

# Create reference prediction data
analysis_tbl <- analysis_tbl %>% 
  mutate(
    pred_df = map(data, ~ tibble(
      age_baseline <-  seq(min(.x$age_baseline, na.rm = TRUE), 
                         max(.x$age_baseline, na.rm = TRUE),
                         length.out = 200),
      sex <-  0,
      edu <-  0,
      gene_apoe <-  factor("e33", levels = levels(.x$gene_apoe)),
      smok_ever <-  0,
    ))
  )

# Produce plots
analysis_tbl <- analysis_tbl %>%
  mutate(
    plot = pmap(list(type, fit, data, pred_df),
               ~ plot_spline(..1, ..2, ..3, ..4, var_name = "age_baseline"))
  )


test_fit <- analysis_tbl$fit[[1]]
test_pred <- tibble(analysis_tbl$pred_df[[1]])

predict(test_fit, newdata = test_pred, type = "lp", se.fit = TRUE)

class(df_ad)

