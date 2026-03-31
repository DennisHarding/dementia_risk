cox_selected <- list("AD" = c("age_baseline", "sex"),
                     "VD" = c("alldm", "is_ih"),
                     "VRD" = c("alldm", "is_ih", "ht"))
legends <- list(
  "gene_apoe_level" = "ApoE Genotype",
  "prs_level" = "Polygenetic Risk Score",
  "sex" = "Sex", 
  "alldm" = "Diabetes", 
  "smok_ever" = "Smoking", 
  "edu" = "Education", 
  "alc" = "Alcohol intake", 
  "is_ih" = "is_ih Stroke",
  "ht" = "Hypertension", 
  "dbp" = "Diastolic BP 10",
  "sbp" = "Systolic BP 10",
  "ihd" = "Ischemic Heart Disease",
  "prs_fac" = "Polygenetic Risk Score",
  "age_baseline" = "Age"
  )
#plot_lasso <- function(type, data){
data <- df_ad
type = "AD"

  X <- data %>%
    select(c("age_baseline", 
             "prs_level", 
             "sex", 
             "gene_apoe_level", 
             "alldm", 
             "smok_ever", 
             "ht", 
             "edu", 
             "alc", 
             "is_ih", 
             "ihd", 
             "dbp10", 
             "sbp10"))
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
    filter((str_detect(variable, str_c(cox_selected[[type]], collapse = "|"))) | 
             (str_detect(variable, "^gene")|str_detect(variable, "^prs")))
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
    filter(!(variable %in% unique(long_coef_df$variable)))

  png(paste0("figures/lasso/lasso_gen_", type, ".png"), 
      width=12, height=6, 
      units="in", 
      res=300)
  par(mar = c(5, 4, 4, 15))
  long_coef_df <- as.data.frame(long_coef_df)
  vars <- unique(long_coef_df$variable)
  lbd <- unique(long_coef_df$log_lambda)
  start_vals <- sapply(vars, function(v) {
    subdf <- long_coef_df[long_coef_df$variable == v, ]
    tail(subdf$coefficient, 1)
  })
  ord <- order(start_vals, decreasing = TRUE)
  vars <- vars[ord]
  plot(NA, 
       xlim = range(long_coef_df$log_lambda), 
       ylim = range(long_coef_df$coefficient),
       xlab = expression(log(lambda)), 
       ylab = "Coefficient")
  
  for (i in seq_along(vars)) {
    subdf <- long_coef_df[long_coef_df$variable == vars[i], ]
    lines(subdf$log_lambda, subdf$coefficient, col = rainbow(length(vars))[i])
  }
  vars
  legend_plot <- unname(unlist(legends[gsub("[0-9]", "", vars)]))
  legend("topright", 
         inset = c(-0.35, 0.3), 
         xpd = TRUE, 
         legend = legend_plot, 
         col = rainbow(length(vars))[seq_along(vars)], 
         lty = 1)
  dev.off()
  
  
  
  png(paste0("figures/lasso/lasso_", type, ".png"), width=12, height=6, units="in", res=300)
  par(mar = c(5, 4, 4, 15))
  long_wo_gen <- as.data.frame(long_wo_gen)
  vars <- unique(long_wo_gen$variable)
  lbd <- unique(long_wo_gen$log_lambda)
  start_vals <- sapply(vars, function(v) {
    subdf <- long_wo_gen[long_wo_gen$variable == v, ]
    tail(subdf$coefficient, 1)
  })
  ord <- order(start_vals, decreasing = TRUE)
  vars <- vars[ord]
  plot(NA, 
       xlim = range(long_wo_gen$log_lambda), 
       ylim = range(long_wo_gen$coefficient),
       xlab = expression(log(lambda)), 
       ylab = "Coefficient")
  
  for (i in seq_along(vars)) {
    subdf <- long_wo_gen[long_wo_gen$variable == vars[i], ]
    lines(subdf$log_lambda, subdf$coefficient, col = rainbow(length(vars))[i])
  }
  legend_plot <-  unlist(ifelse(gsub("[0-9]+$", "", vars) %in% names(legends), legends[gsub("[0-9]+$", "", vars)], gsub("[0-9]+$", "", vars)))
  legend("topright", 
         inset = c(-0.35, 0.3), 
         xpd = TRUE, 
         legend = legend_plot, 
         col = rainbow(length(vars))[seq_along(vars)], 
         lty = 1)
  dev.off()
  
  cvfit <- cv.glmnet(X_clean, y_clean, family = "cox")
  
  plot(cvfit)
  
}

