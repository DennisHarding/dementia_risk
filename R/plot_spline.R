plot_spline <- function(demtype, fit, data, pred_df, var_name) {
  
  var_name = as.character(var_name)
  
  # Predict
  pred <- predict(fit, newdata = pred_df, type = "lp", se.fit = TRUE)
  
  pred_df$hr <- exp(pred$fit - mean(pred$fit))
  pred_df$lower <- exp(pred$fit - 1.96 * pred$se.fit - mean(pred$fit))
  pred_df$upper <- exp(pred$fit + 1.96 * pred$se.fit - mean(pred$fit))
  
  # Density of chosen variable
  dens <- density(data[[var_name]], na.rm = TRUE)
  dens_df <- data.frame(x = dens$x, y = dens$y)
  
  hr_range <- range(pred_df$lower, pred_df$upper)
  dens_scale <- 0.2 * (hr_range[2] - hr_range[1]) / max(dens_df$y)
  dens_df$y_scaled <- hr_range[1] + dens_scale * dens_df$y

  # Plot
  ggplot() +
    geom_line(data = pred_df, aes(x = .data[[var_name]], y = hr)) +
    geom_ribbon(data = pred_df,
                aes(x = .data[[var_name]], ymin = lower, ymax = upper),
                alpha = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_line(data = dens_df,
              aes(x = x, y = y_scaled),
              linewidth = 0.3) +
    geom_ribbon(data = dens_df,
                aes(x = x, ymin = hr_range[1], ymax = y_scaled),
                fill = "steelblue", alpha = 0.3) +
    labs(x = var_name,
         y = "Hazard Ratio",
         title = paste("HR as a function of", var_name,"for", demtype)) +
    theme_minimal()
}
