fit_cox_spline <- function(data, vars, var, df = 3) {
  
  if (!var %in% vars) {
    stop("var must be included in vars")
  }
  vars <- setdiff(vars, var)

  spline_term <- paste0("ns(", var, ", df = ", df,")")
  
  rhs <- paste(c(vars, spline_term), collapse = " + ")

  formula <- as.formula(
    paste0("Surv(futime, fail_bin) ~ ", rhs)
  )
  
  fit <- survival::coxph(formula, data = data)
  return(fit)
}



create_pred_df <- function(data, var) {
  
  var_seq <- seq(
    min(data[[var]], na.rm = TRUE),
    max(data[[var]], na.rm = TRUE),
    length.out = 200
  )
  
  base_row <- tibble(
    prs = mean(data$prs, na.rm = TRUE),
    sex = factor(0, levels = levels(data$sex)),
    gene_apoe = factor("e33", levels = levels(data$gene_apoe)),
    smok_ever = factor(0, levels = levels(data$smok_ever)),
    alc = factor(0, levels = levels(data$alc)),
    ht = factor(0, levels = levels(data$ht)),
    is_ih = factor(0, levels = levels(data$is_ih)),
    ihd = factor(0, levels = levels(data$ihd)),
    edu_cont = mean(data$edu_cont, na.rm = TRUE),
    age_baseline = mean(data$age_baseline, na.rm = TRUE),
    dbp10 = mean(data$dbp10, na.rm = TRUE),
    sbp10 = mean(data$sbp10, na.rm = TRUE)
  )
  
  pred_df <- base_row %>%
    slice(rep(1, 200)) %>%
    mutate(!!var := var_seq)
  
  return(pred_df)
}

plot_spline_dfs <- function(type, data, vars, var, dfmin, dfmax, save = FALSE) {
  
  df_seq <- seq(dfmin, dfmax)
  
  fits <- map(df_seq, function(d) {
    fit_cox_spline(
      data = data,
      vars = vars,
      var = var,
      df = d
    )
  })
  
  names(fits) <- paste0("df", df_seq)
  
  pred_df <- create_pred_df(data, var)
  
  var = as.character(var)
  
  ref_model <- fits[[1]]   # smallest df
  ref_pred <- predict(ref_model, newdata = pred_df, type = "lp")
  ref_val <- mean(ref_pred)
  
  
  pred_dfs <- imap(fits, function(fit, name) {
    pred <- predict(fit, newdata = pred_df, type = "lp", se.fit = TRUE)
    df_out <- pred_df
    df_out$hr <- pred$fit - ref_val
    df_out$lower <- pred$fit - 1.96*pred$se.fit - ref_val
    df_out$upper <- pred$fit + 1.96*pred$se.fit - ref_val
    df_out$model <- name
    df_out
  })
  pred_all <- bind_rows(pred_dfs)
  
  # Density of chosen variable
  dens <- density(data[[var]], na.rm = TRUE)
  dens_df <- data.frame(x = dens$x, y = dens$y)
  
  hr_range <- range(pred_all$lower, pred_all$upper)
  dens_scale <- 0.2 * (hr_range[2] - hr_range[1]) / max(dens_df$y)
  dens_df$y_scaled <- hr_range[1] + dens_scale * dens_df$y
  
  # Plot
  p <- ggplot(pred_all, aes(x = .data[[var]], y = hr, color = model)) +
    geom_ribbon(linetype = "dashed", linewidth = 0.25, aes(ymin = lower, ymax = upper),
                alpha = 0.05) +
    geom_line(linewidth = 0.7) +
    geom_line(data = dens_df, 
              aes(x = x, y = y_scaled), 
              inherit.aes = FALSE, linewidth = 0.3) +
    geom_ribbon(data = dens_df,
              aes(x = x, ymin = hr_range[1], ymax = y_scaled),
              inherit.aes =  FALSE,
              fill = "steelblue", alpha = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(
      x = var,
      y = "log Hazard Ratio",
      title = paste("log HR as a function of", var, "for", type)) +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  if (save == TRUE){
    
    filename = glue("multisplineplot_{type}_{var}_{dfmin}{dfmax}.png")
    path = "figures/splineplot/"
    
    print(glue("Saving {filename} into {path}"))
    
    ggsave(
      filename,
      width = 5, 
      height = 5, 
      plot = p,
      path = path
    )
  }
  
  return(p)
}


