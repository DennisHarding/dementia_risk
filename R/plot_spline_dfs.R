
create_pred_df <- funciton(data, var){
  pred_df <- tibble(
    prs = mean(data$prs, na.rm = TRUE),
    sex = 0,
    edu = 0,
    gene_apoe = factor("e33", levels = levels(data$gene_apoe)),
    smok_ever = 0,
    alc = 0,
    stroke = 0,
    
    edu_cont = mean(data$edu_cont, na.rm = TRUE),
    age_baseline = mean(data$age_baseline, na.rm = TRUE),
    dbp10 = mean(data$dbp10, na.rm = TRUE),
    sbp10 = mean(data$sbp10, na.rm = TRUE)
  )
  pred_df %>% mutate(
    var = seq(min(data[[var]], na.rm = TRUE),
          max(data[[var]], na.rm = TRUE),
          length.out = 200)

  )
}


analysis_tbl <- analysis_tbl %>% 
  mutate(
    pred_df = map(data, ~ tibble(
      prs = mean(.x$prs, na.rm = TRUE),
      sex = 0,
      edu = 0,
      gene_apoe = factor("e33", levels = levels(.x$gene_apoe)),
      smok_ever = 0,
      alc = 0,
      stroke = 0,
    
      edu_cont = mean(.x$edu_cont, na.rm = TRUE),
      age_baseline = mean(.x$age_baseline, na.rm = TRUE),
      dbp10 = mean(.x$dbp10, na.rm = TRUE),
      sbp10 = mean(.x$sbp10, na.rm = TRUE)
    ))
  )


var = age_baseline

start <- 1
stop <- 4

vars <- paste0("df", start:stop)
vars
spline_tbl <- setNames(as.list(vars), start:stop)

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
        smok_ever +
        alc +
        edu_cont +
        dbp10 +
        sbp10
      ,
      data = .x
    ))
  )
# Generate reference prediction data

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