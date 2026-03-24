#> -----------------------------
#> LOAD LIBRARIES
#> -----------------------------
{
library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(glue)
library(gtable)
library(patchwork)
library(broom)
library(glmnet)
library(forestplot)
}
#> -----------------------------
#> DATA FORMATTING
#> -----------------------------
{data <- fread("data/custom_data_1.tsv")

source("R/formatting.R")

df <- format_df(data)

df_ad <- format_df_ad(df)
df_vd <- format_df_vd(df)
df_vrd <- format_df_vrd(df)
}
#> -----------------------------
#> Prepare analysis table
#> -----------------------------
{
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
}
#> -----------------------------
#> Forest plots
#> -----------------------------
{
source("R/plot_forest.R")

analysis_tbl = analysis_tbl %>% 
  mutate(forest_plot = pmap(list(type, data), ~ 
                              plot_forest(..1, ..2)))
}
#> -----------------------------
#> DEFINE & FIT COX MODELS
#> -----------------------------
{analysis_tbl <- analysis_tbl %>%
  mutate(
    data = map(data, ~ filter(.x, !is.na(sbp),
                              !is.na(dbp),
                              !is.na(ht))
    )
  )
analysis_tbl <- analysis_tbl %>% 
  mutate(
    cox_fit = map(data, ~ coxph(
      Surv(futime, fail_bin) ~ 
        age_baseline +
        sex + 
        prs +
        gene_apoe +
        edu_cont +
        smok_ever +
        alc +
        dbp10 +
        sbp10
      ,
      data = .x
    ))
  )

analysis_tbl %>% 
  mutate(cox_summary = map(cox_fit, summary)) %>% 
  pull(cox_summary) %>% 
  walk(print)
# Optional: define candidate cox model:

analysis_tbl <- analysis_tbl %>% 
  mutate(
    cox_fit_new = map(data, ~ coxph(
      Surv(futime, fail_bin) ~ 
        I(age_baseline^2) +
        sex + 
        prs +
        gene_apoe +
        edu_cont +
        smok_ever +
        alc +
        ht

      ,
      data = .x
    ))
  )


analysis_tbl %>% 
  mutate(cox_summary_new = map(cox_fit_new, summary)) %>% 
  pull(cox_summary_new) %>% 
  walk(print)
}
#> -----------------------------
#> Anova test
#> -----------------------------
{
analysis_tbl <- analysis_tbl %>%
  mutate(
    anova_test = pmap(list(cox_fit, cox_fit_new), ~ 
      anova(..1,..2)
    )
  )

# Print anova tests
analysis_tbl %>% pull(anova_test) %>% walk(print)
}
#> -----------------------------
#> AIC test
#> -----------------------------
{
analysis_tbl <- analysis_tbl %>%
  mutate(
    AIC_test = pmap(list(cox_fit, cox_fit_new), ~ 
      AIC(..1,..2)
    )
  )

# Print AIC tests
analysis_tbl %>% pull(AIC_test) %>% walk(print)

}
#> -----------------------------
#> SPLINES
#> -----------------------------
{source("R/plot_spline.R")

analysis_tbl <- analysis_tbl %>% 
  mutate(
    cox_fit_spline = map(data, ~ coxph(
      Surv(futime, fail_bin) ~ 
        age_baseline +
        sex + 
        prs +
        gene_apoe +
        edu_cont +
        smok_ever +
        alc +
        dbp10 +
        sbp10
      ,
      data = .x
    ))
  )
# Generate reference prediction data
analysis_tbl <- analysis_tbl %>% 
  mutate(
    pred_df = map(data, ~ tibble(
      prs = mean(.x$prs, na.rm = TRUE),
      sex = 0,
      edu = 0,
      gene_apoe = factor("e33", levels = levels(.x$gene_apoe)),
      smok_ever = 0,
      alc = 0,
      edu_cont = mean(.x$edu_cont, na.rm = TRUE),
      stroke = 0,
      age_baseline = mean(.x$age_baseline, na.rm = TRUE),
      dbp = 0,
      sbp = seq(min(.x$sbp, na.rm = TRUE),
                max(.x$sbp, na.rm = TRUE),
                length.out = 200)
    ))
  )

# Produce spline plots
analysis_tbl <- analysis_tbl %>%
  mutate(
    spline_plot = pmap(list(type, cox_fit_spline, data, pred_df),
               ~ plot_spline(..1, ..2, ..3, ..4, 
                             var = "sbp",
                             save = TRUE))
  )

# PRINT SPLINE PLOTS
analysis_tbl %>% pull(spline_plot) %>% walk(print)

# # Produce spline plots WITH COMPARISON
# analysis_tbl <- analysis_tbl %>%
#   mutate(
#     spline_plot_vs = pmap(list(type, cox_fit, data, pred_df, cox_fit_new),
#                ~ plot_spline(..1, ..2, ..3, ..4, ..5, 
#                              var_name = "age_baseline",
#                              save = FALSE)
#     )
#   )
# 
# 
# # Produce spline plots
# analysis_tbl %>% pull(spline_plot_vs) %>% walk(print)
}
#> -----------------------------
#> Proportional Hazards tests - Schoenfeld residuals
#> -----------------------------
{source("R/plot_sfr.R")

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
    ph_plot = pmap(list(type, ph_test), ~ 
      plot_sfr(..1, ..2, save = FALSE))
  )

# PRINT PH PLOTS
analysis_tbl %>% pull(ph_plot) %>% walk(print)
}
#> -----------------------------
#> Proportional Hazards tests - Log-Log plots
#> -----------------------------
{source("R/plot_loglog.R")

analysis_tbl <- analysis_tbl %>%
  mutate(surv_curvs = pmap(list(type, data), ~ {
    
    # SPECIFY vars of interest -->
    vars = c("sbp", "dbp")
    
    plot_loglog(..1, ..2, vars, save = FALSE)
  }))

analysis_tbl %>% pull(surv_curvs) %>% walk(print)
}
#> -----------------------------
#> Generate Summary Table
#> -----------------------------
{
library(gtsummary)

#remember to select variables to include in summary table

summary_table <- 
  df %>%
  select(
    "age_baseline",
    "sex",
    "edu_cont",
    "edu",
    "gene_apoe",
    "alldem_bin",
    "ad_bin",
    "vd_bin",
    "vrd_bin",
    "prs",
    "alldm",
    "sbp",
    "dbp"
  ) %>%
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
}