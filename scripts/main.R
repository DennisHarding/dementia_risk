library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(glue)

#> -----------------------------
#> DATA FORMATTING
#> -----------------------------
source("R/formatting.R")

data <- fread("data/custom_data_1.tsv")

df <- format_df(data)

df_ad <- format_df_ad(df)
df_vd <- format_df_vd(df)
df_vrd <- format_df_vrd(df)
#df_alldem <- format_df_alldem()

#> -----------------------------
#> Prepare analysis table
#> -----------------------------

types_list <- list(
  AD = "AD",
  VD = "VD",
  VRD = "VRD"
)

data_list <- list(
  AD = df_ad,
  VD = df_vd,
  VRD = df_vrd
)

analysis_tbl <- tibble(
  type = types_list,
  data = data_list
)


#> prs:gene_apoe34 significant interaction for VD
#> 

#> -----------------------------
#> DEFINE & FIT COX MODELS
#> -----------------------------

#> I want to be able to compare cox models with different variables,
#> issue is that different variables have differnet missing values
#> and i dont want to simply remove all observations which have missing
#> values,
#> 
#> instead i want to remove only the values which are in the columns 
#> which are different between the two cox models. So i could make a 
#> simple row of code in which the user can take the data stored in 
#> analysis_tbl and filter the datasets to only keep obsevations which
#> are not na in the differing variabble.

analysis_tbl <- analysis_tbl %>%
  mutate(
    data = map(data, ~ filter(.x, !is.na(hf))
    )
  )
analysis_tbl
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
        ht +
        hf
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
source("R/plot_spline.R")

# Generate reference prediction data
analysis_tbl <- analysis_tbl %>% 
  mutate(
    pred_df = map(data, ~ tibble(
      prs = mean(.x$prs, na.rm = TRUE),
      sex = 0,
      edu = 0,
      gene_apoe = factor("e33", levels = levels(.x$gene_apoe)),
      smok_ever = 0,
      age_baseline = seq(min(.x$age_baseline, na.rm = TRUE),
                max(.x$age_baseline, na.rm = TRUE),
                length.out = 200),
      ht = 0,
      hf = 0
    ))
  )

# Produce spline plots
analysis_tbl <- analysis_tbl %>%
  mutate(
    spline_plot = pmap(list(type, cox_fit, data, pred_df),
               ~ plot_spline(..1, ..2, ..3, ..4, 
                             var_name = "age_baseline"))
  )

# Produce spline plots WITH COMPARISON
analysis_tbl <- analysis_tbl %>%
  mutate(
    spline_plot = pmap(list(type, cox_fit, data, pred_df, cox_fit_new),
               ~ plot_spline(..1, ..2, ..3, ..4, ..5, 
                             var_name = "age_baseline"))
  )

# PRINT SPLINE PLOTS
analysis_tbl %>% pull(spline_plot) %>% walk(print)


analysis_tbl
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
library(ggforestplot)
source("R/plot_forest.R")

analysis_tbl <- analysis_tbl %>%
  mutate(
    forest_plot = pmap(list(type, cox_fit), ~
      plot_forest(..1, ..2))
  )

# PRINT FOREST PLOTS
analysis_tbl %>% pull(forest_plot) %>% walk(print)

#> -----------------------------
#> Generate Summary Table
#> -----------------------------

library(gtsummary)

#select variables to include in summary table
vars <- c("age_baseline",
          "sex",
          "edu_yrs",
          "gene_apoe",
          "alldem_bin",
          "prs",
          "ad_bin",
          "vd_bin",
          "vrd_bin",
          "alldm"
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
