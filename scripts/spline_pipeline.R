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

# Required scripts and functions
source("scripts/formatting.R")
source("R/plot_spline.R")

# Reference variables and data frames > analysis table
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

# COX MODELS
analysis_tbl <- analysis_tbl %>% 
  mutate(
    cox_fit = map(data, ~ coxph(
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

# Reference prediction data
analysis_tbl <- analysis_tbl %>% 
  mutate(
    pred_df = map(data, ~ tibble(
      age_baseline = seq(min(.x$age_baseline, na.rm = TRUE), 
                         max(.x$age_baseline, na.rm = TRUE),
                         length.out = 200),
      sex = 0,
      edu = 0,
      gene_apoe = factor("e33", levels = levels(.x$gene_apoe)),
      smok_ever = 0,
    ))
  )

# Produce spline plots
analysis_tbl <- analysis_tbl %>%
  mutate(
    spline_plot = pmap(list(type, cox_fit, data, pred_df),
               ~ plot_spline(..1, ..2, ..3, ..4, var_name = "age_baseline"))
  )

# Schoenfeld residual tests (proportional hazards tests)
analysis_tbl <- analysis_tbl %>%
  mutate(
    ph_test = map(cox_fit, ~
      cox.zph(.x)
    )
  )

# Plotting Schoenfeld residuals
analysis_tbl <- analysis_tbl %>%
  mutate(
    ph_plot = map(ph_test, ~
      ggcoxzph(.x, ggtheme = theme_bw())
    )
  )

# PRINT SPLINE PLOTS
analysis_tbl %>% pull(spline_plot) %>% walk(print)

# PRINT PH PLOTS
analysis_tbl %>% pull(ph_plot) %>% walk(print)
