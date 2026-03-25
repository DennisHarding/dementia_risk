cox_selected <- list("AD" = c("edu", "ht"),
                     "VD" = c("alldm", "is_ih"),
                     "VRD" = c("alldm", "is_ih", "ht"))
legends <- list(
  "sex" = "Sex", 
  "alldm" = "Diabetes", 
  "smok_ever" = "Smoking", 
  "edu" = "Education", 
  "alc" = "Alcohol intake", 
  "is_ih" = "is_ih Stroke",
  "ht" = "Hypertension", 
  "dbp10" = "Diastolic BP 10",
  "sbp10" = "Systolic BP 10",
  "ihd" = "Ischemic Heart Disease",
  "prs_fac" = "Polygenetic Risk Score"
  
  
  )
# 
# df <- df %>%
#   select(t2dm_bin, smok_ever, iscvd_bin, iscvd_date_first, giga_ih_bin, giga_ih_date_first, ht_bin, ht_date_first, edu_grp, alc_freq, pa_met_vigorous, pa_met_moderate, gene_apoe, ad_bin, nonad_bin, vd_bin, alldem_bin, phylo_score, birth_date, sex, ad_date_first, nonad_date_first, vd_date_first, alldem_date_first, t2dm_date_first, entry_date, chol, hdl) %>%
#   mutate(
#     age_baseline = as.numeric(difftime(entry_date, birth_date, unit = "days")) / 365.25,
#     age_ad_bin = as.numeric(difftime(ad_date_first, birth_date, unit = "days")) / 365.25,
#     age_nonad_bin = as.numeric(difftime(nonad_date_first, birth_date, unit = "days")) / 365.25, 
#     age_vd_bin = as.numeric(difftime(vd_date_first, birth_date, unit = "days")) / 365.25,
#     age_alldem_bin = as.numeric(difftime(alldem_date_first, birth_date, unit = "days")) / 365.25,
#     polygenic = ntile(phylo_score, 4),
#     education = case_when(str_detect(edu_grp, "1") ~ 1,
#                           str_detect(edu_grp, "5") ~ 1,
#                           str_detect(edu_grp,"6") ~ 1,
#                           str_detect(edu_grp,"2") ~ 1,
#                           str_detect(edu_grp, "-3") ~ NA, 
#                           str_detect(edu_grp, "3|4") ~ 0,
#                           str_detect(edu_grp, "-7") ~ 0,
#                           .default = NA
#     ),
#     smok_ever = case_when(
#       smok_ever == "Yes" ~ 1,
#       smok_ever == "No" ~ 0,
#       .default = NA
#     ),
#     alc_freq = case_when(
#       alc_freq %in% c(1,2) ~ 1,
#       alc_freq %in% c(3,4,5,6) ~ 0,
#       .default = NA
#     ),
#     physical = case_when(
#       pa_met_vigorous >= 75 | pa_met_moderate >= 150 ~ 1,
#       pa_met_vigorous < 75 & pa_met_moderate < 150 ~ 0,
#       .default = NA
#     ),
#     nonhdl = ntile(chol - hdl, 2),
#     stroke = case_when(
#       (iscvd_bin == 1 & iscvd_date_first <= entry_date) | (giga_ih_bin == 1 & giga_ih_date_first <= entry_date) ~ 1,
#       (iscvd_bin == 1 & iscvd_date_first > entry_date) & (giga_ih_bin == 1 & giga_ih_date_first > entry_date) ~ 0,
#       (iscvd_bin == 1 & iscvd_date_first > entry_date) & giga_ih_bin == 0 ~ 0,
#       iscvd_bin == 0 & (giga_ih_bin == 1 & giga_ih_date_first > entry_date) ~ 0,
#       iscvd_bin == 0 & giga_ih_bin == 0 ~ 0,
#       .default = NA
#     ),
#     ht = case_when(
#       (ht_bin == 1) & (ht_date_first <= entry_date) ~ 1,
#       (ht_bin == 1) & (ht_date_first > entry_date) ~ 0,
#       ht_bin == 0 ~ 0,
#       .default = NA
#     ),
#     t2dm = case_when(
#       (t2dm_bin == 1) & (t2dm_date_first <= entry_date) ~ 1,
#       (t2dm_bin == 1) & (t2dm_date_first > entry_date) ~ 0,
#       t2dm_bin == 0 ~ 0,
#       .default = NA
#     )
#   )

plot_lasso <- function(type, data){
  X <- data %>%
    select(c("age_baseline", "prs_fac", "sex", "gene_apoe", "alldm", "smok_ever", "ht", "edu", "alc", "is_ih", "ihd", "dbp10", "sbp10"))
  y <- Surv(data$futime, event = data$fail_bin)
  complete_cases <- complete.cases(X, y)
  X_clean <- model.matrix(~ ., data = X[complete_cases, ])[,-1]
  y_clean <- y[complete_cases]
  fit <- glmnet(x = X_clean, y = y_clean, family = "cox", alpha = 1) 
  
  lambda_seq <- fit$lambda
  coefs <- matrix(NA, nrow = nrow(coef(fit)), ncol = length(lambda_seq))
  rownames(coefs) <- rownames(coef(fit))
  for (i in seq_along(lambda_seq)) {
    coefs[, i] <- as.vector(coef(fit, s = lambda_seq[i]))
  }
  
  rownames(coefs) <- rownames(coef(fit))
  colnames(coefs) <- paste0("lambda_", seq_along(lambda_seq))
  coef_df <- as.data.frame(coefs)
  coef_df$variable <- rownames(coef_df)
  
  long_coef_df <- coef_df %>%
    pivot_longer(
      cols = starts_with("lambda_"),
      names_to = "lambda_id",
      values_to = "coefficient"
    ) %>%
    mutate(
      lambda_index = as.integer(gsub("lambda_", "", lambda_id)),
      lambda = lambda_seq[lambda_index],
      log_lambda = log(lambda)
    ) %>%
    filter((str_detect(variable, str_c(cox_selected[[type]], collapse = "|"))) | (str_detect(variable, "^gene")|str_detect(variable, "^prs")))
  long_wo_gen <- coef_df %>%
    pivot_longer(
      cols = starts_with("lambda_"),
      names_to = "lambda_id",
      values_to = "coefficient"
    ) %>%
    mutate(
      lambda_index = as.integer(gsub("lambda_", "", lambda_id)),
      lambda = lambda_seq[lambda_index],
      log_lambda = log(lambda)
    ) %>%
    filter(!(str_detect(variable, "^gene")) & !(str_detect(variable, "^prs_fac")))
  
  png(paste0("figures/lasso/lasso_gen_", type, ".png"), width=12, height=6, units="in", res=300)
  par(mar = c(5, 4, 4, 12))
  long_coef_df <- as.data.frame(long_coef_df)
  vars <- unique(long_coef_df$variable)
  lbd <- unique(long_coef_df$log_lambda)
  start_vals <- sapply(vars, function(v) {
    subdf <- long_coef_df[long_coef_df$variable == v, ]
    tail(subdf$coefficient, 1)
  })
  ord <- order(start_vals, decreasing = TRUE)
  vars <- vars[ord]
  plot(NA, xlim = range(long_coef_df$log_lambda), ylim = range(long_coef_df$coefficient),
       xlab = expression(log(lambda)), ylab = "Coefficient")
  for (i in seq_along(vars)) {
    subdf <- long_coef_df[long_coef_df$variable == vars[i], ]
    lines(subdf$log_lambda, subdf$coefficient, col = rainbow(length(vars))[i])
  }
  legend_plot <-  unlist(ifelse(gsub("[0-9]+$", "", vars) %in% names(legends), legends[gsub("[0-9]+$", "", vars)], gsub("[0-9]+$", "", vars)))
  legend("topright", inset = c(-0.25, 0.3), xpd = TRUE, legend = legend_plot, col = rainbow(length(vars))[seq_along(vars)], lty = 1)
  dev.off()
  
  png(paste0("figures/lasso/lasso_", type, ".png"), width=12, height=6, units="in", res=300)
  par(mar = c(5, 4, 4, 12))
  long_wo_gen <- as.data.frame(long_wo_gen)
  vars <- unique(long_wo_gen$variable)
  lbd <- unique(long_wo_gen$log_lambda)
  start_vals <- sapply(vars, function(v) {
    subdf <- long_wo_gen[long_wo_gen$variable == v, ]
    tail(subdf$coefficient, 1)
  })
  ord <- order(start_vals, decreasing = TRUE)
  vars <- vars[ord]
  plot(NA, xlim = range(long_wo_gen$log_lambda), ylim = range(long_wo_gen$coefficient),
       xlab = expression(log(lambda)), ylab = "Coefficient")
  for (i in seq_along(vars)) {
    subdf <- long_wo_gen[long_wo_gen$variable == vars[i], ]
    lines(subdf$log_lambda, subdf$coefficient, col = rainbow(length(vars))[i])
  }
  legend_plot <-  unlist(ifelse(gsub("[0-9]+$", "", vars) %in% names(legends), legends[gsub("[0-9]+$", "", vars)], gsub("[0-9]+$", "", vars)))
  legend("topright", inset = c(-0.25, 0.3), xpd = TRUE, legend = legend_plot, col = rainbow(length(vars))[seq_along(vars)], lty = 1)
  dev.off()
}