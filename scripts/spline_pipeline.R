install.packages("survival")
library(survival)
install.packages("survminer")
library(survminer)
library(tibble)
library(dplyr)
library(purrr)

source("splines/plot_spline.R")

ph_test <- cox.zph(fit)

ggcoxzph(ph_test, ggtheme = theme_bw())



data_list <- list(
  AD = df_ad,
  VD = df_vd,
  NonAD = df_nonad
)

analysis_tbl <- tibble(
  type = names(data_list),
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
    pred_df = map(data, ~ {
      age_baseline <-  seq(min(df$age_baseline, na.rm = TRUE), 
                         max(df_ad$age_baseline, na.rm = TRUE),
                         length.out = 200)
      sex <-  0
      edu <-  0
      gene_apoe <-  factor("e33", levels = levels(df_ad$gene_apoe))
      smok_ever <-  0
    })
  )

# Produce plots
analysis_tbl <- analysis_tbl %>%
  mutate(
    plot = pmap(list(fit, data, pred_df),
               ~ plot_spline(..1, ..2, ..3, var_name = "age_baseline"))
  )



# Fit cox model
fit <- coxph(Surv(futime_ad, ad_bin) ~ 
               age_baseline +
               sex + 
               edu +
               gene_apoe +
               smok_ever # + 
              ,
             data = df_ad)


pred_df_ad <- tibble(
  age_baseline = seq(min(df_ad$age_baseline), 
            max(df_ad$age_baseline),
            length.out = 200),
  sex = 0,
  edu = 0,
  gene_apoe = factor("e33", levels = levels(df_ad$gene_apoe)),
  smok_ever = 0
)
pred_df_vd <- tibble(
  age_baseline = seq(min(df_vd$age_baseline), 
            max(df_vd$age_baseline),
            length.out = 200),
  sex = 0,
  edu = 0,
  gene_apoe = factor("e33", levels = levels(df_vd$gene_apoe)),
  smok_ever = 0
)

pred_df_nonad <- tibble(
  age_baseline = seq(min(df_nonad$age_baseline), 
            max(df_nonad$age_baseline),
            length.out = 200),
  sex = 0,
  edu = 0,
  gene_apoe = factor("e33", levels = levels(df_nonad$gene_apoe)),
  smok_ever = 0
)


plot_cox_curve(fit, df_ad, pred_df, "age_baseline")
plot_cox_curve(fit, df_vd, pred_df, "age_baseline")
plot_cox_curve(fit, df_nonad, pred_df, "age_baseline")

