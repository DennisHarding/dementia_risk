install.packages("survival")
library(survival)
install.packages("survminer")
library(survminer)
library(splines)
source("splines/plot_spline.R")

fit <- coxph(
  Surv(futime, fail_bin) ~ 
    gene_apoe +
    prs +
    alldm +
    ns(age_baseline,df = 1) +
    sex +
    edu +
    smok_ever +
    ht +
    non_hdl
    ,
  data = df_alldem)
ph_test <- cox.zph(fit)

ggcoxzph(ph_test, ggtheme = theme_bw())


pred_df <- tibble(
  gene_apoe = factor("e33", levels = levels(df_alldem$gene_apoe)),
  prs = median(df_alldem$prs, na.rm = TRUE),
  alldm = 0,
  age_baseline = median(df_alldem$age_baseline),
  sex = 0,
  edu = 0,
  smok_ever = 0,
  ht = 0,
  non_hdl = seq(min(df_alldem$non_hdl, na.rm = TRUE),
                max(df_alldem$non_hdl, na.rm = TRUE),
                length.out = 200)
)


plot_spline(fit, df_alldem, pred_df, "non_hdl")