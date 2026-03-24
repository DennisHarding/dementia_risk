
paste0(df, 1:5)


spline_tbl = list(
  df1 = "df1",
  df2 = "df2",
  df3 = "df3",
  df4 = "df4"
)


names = c(
  "age_baseline",
  "sex",
  "prs",
  "gene_apoe",
  "edu_cont",
  "smok_ever",
  "alc",
  "dbp10",
  "sbp10")

spline_tbl <- spline_tbl %>% 
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