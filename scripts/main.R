library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(ggforestplot)
library(glue)

# Required scripts and functions
source("R/formatting.R")
source("R/plot_spline.R")
source("R/plot_forest.R")

#> -----------------------------
#> DATA FORMATTING
#> -----------------------------

data <- fread("data/custom_data_1.tsv")

df <- format_df(data)

df_ad <- format_df_ad(df)
df_vd <- format_df_vd(df)
df_nonad <- format_df_nonad(df)
#df_alldem <- format_df_alldem()

#> -----------------------------
#> Prepare analysis table
#> -----------------------------
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

#> -----------------------------
#> DEFINE & FIT COX MODELS
#> -----------------------------
analysis_tbl <- analysis_tbl %>% 
  mutate(
    cox_fit = map(data, ~ coxph(
      Surv(futime, fail_bin) ~ 
        I(age_baseline^2) +
        sex + 
        edu +
        smok_ever +
        prs +
        gene_apoe +
        ht
      ,
      data = .x
    ))
  )

# Optional: define candidate cox model:

analysis_tbl <- analysis_tbl %>% 
  mutate(
    cox_fit_new = map(data, ~ coxph(
      Surv(futime, fail_bin) ~ 
        I(age_baseline^2) +
        sex + 
        edu +
        smok_ever +
        prs +
        gene_apoe +
        ht
      ,
      data = .x
    ))
  )


# Print summary of cox models

analysis_tbl %>% 
  mutate(cox_summary = map(cox_fit, summary)) %>% 
  pull(cox_summary) %>% 
  walk(print)

analysis_tbl %>% 
  mutate(cox_summary_new = map(cox_fit_new, summary)) %>% 
  pull(cox_summary_new) %>% 
  walk(print)

#> -----------------------------
#> Anova test
#> -----------------------------
analysis_tbl <- analysis_tbl %>%
  mutate(
    anova_test = pmap(list(cox_fit, cox_fit_new), ~ 
      anova(..1,..2)
    )
  )

# Print anova tests
analysis_tbl %>% pull(anova_test) %>% walk(print)

#> -----------------------------
#> AIC test
#> -----------------------------
analysis_tbl <- analysis_tbl %>%
  mutate(
    AIC_test = pmap(list(cox_fit, cox_fit_new), ~ 
      AIC(..1,..2)
    )
  )

# Print AIC tests
analysis_tbl %>% pull(AIC_test) %>% walk(print)


#> -----------------------------
#> SPLINES
#> -----------------------------

# Generate reference prediction data
analysis_tbl <- analysis_tbl %>% 
  mutate(
    pred_df = map(data, ~ tibble(
      age_baseline = mean(.x$age_baseline, na.rm = TRUE),
      sex = 0,
      edu = 0,
      gene_apoe = factor("e33", levels = levels(.x$gene_apoe)),
      smok_ever = 0,
      prs = seq(min(.x$prs, na.rm = TRUE),
                max(.x$prs, na.rm = TRUE),
                length.out = 200),
      ht = 0
    ))
  )

# Produce spline plots
analysis_tbl <- analysis_tbl %>%
  mutate(
    spline_plot = pmap(list(type, cox_fit, data, pred_df),
               ~ plot_spline(..1, ..2, ..3, ..4, var_name = "prs"))
  )

# PRINT SPLINE PLOTS
analysis_tbl %>% pull(spline_plot) %>% walk(print)

#> -----------------------------
#> Proportional Hazards tests (Schoenfeld residuals)
#> -----------------------------

# Tests
analysis_tbl <- analysis_tbl %>%
  mutate(
    ph_test = map(cox_fit, ~
      cox.zph(.x)
    )
  )

# Plots
analysis_tbl <- analysis_tbl %>%
  mutate(
    ph_plot = pmap(list(type, ph_test), ~{
      caption = glue("PH Test For  {..1}")
      ggcoxzph(caption = caption, fit = ..2, ggtheme = theme_bw())
    })
  )

# PRINT PH PLOTS
analysis_tbl %>% pull(ph_plot) %>% walk(print)

#> -----------------------------
#> Forest plots
#> -----------------------------
analysis_tbl <- analysis_tbl %>%
  mutate(
    forest_plot = pmap(list(type, cox_fit), ~
      plot_forest(..1, ..2))
  )

# PRINT FOREST PLOTS
analysis_tbl %>% pull(forest_plot) %>% walk(print)

