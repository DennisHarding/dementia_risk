plot_spline <- function(type, fit, data, pred_df, fit_new = NULL, var, save = FALSE) {
  
  var = as.character(var)
  
  # Predict
  pred <- predict(fit, newdata = pred_df, type = "lp", se.fit = TRUE)
  
  pred_df$hr <- exp(pred$fit - mean(pred$fit))
  pred_df$lower <- exp(pred$fit - 1.96 * pred$se.fit - mean(pred$fit))
  pred_df$upper <- exp(pred$fit + 1.96 * pred$se.fit - mean(pred$fit))
  
  # Predict new
  if (!is.null(fit_new)){
    pred_new <- predict(fit_new, newdata = pred_df, type = "lp", se.fit = TRUE)
    
    pred_df$hr_new <- exp(pred_new$fit - mean(pred_new$fit))
    pred_df$lower_new <- exp(pred_new$fit - 1.96 * pred_new$se.fit - mean(pred_new$fit))
    pred_df$upper_new <- exp(pred_new$fit + 1.96 * pred_new$se.fit - mean(pred_new$fit))
  }
  # Density of chosen variable
  dens <- density(data[[var]], na.rm = TRUE)
  dens_df <- data.frame(x = dens$x, y = dens$y)
  
  hr_range <- range(pred_df$lower, pred_df$upper)
  dens_scale <- 0.2 * (hr_range[2] - hr_range[1]) / max(dens_df$y)
  dens_df$y_scaled <- hr_range[1] + dens_scale * dens_df$y

  # Plot
  p <- ggplot() +
    geom_line(data = pred_df, aes(x = .data[[var]], y = hr)) +
    geom_ribbon(data = pred_df,
                aes(x = .data[[var]], ymin = lower, ymax = upper),
                alpha = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_line(data = dens_df,
              aes(x = x, y = y_scaled),
              linewidth = 0.3) +
    geom_ribbon(data = dens_df,
                aes(x = x, ymin = hr_range[1], ymax = y_scaled),
                fill = "steelblue", alpha = 0.3) +
    labs(x = var,
         y = "Hazard Ratio",
         title = paste("HR as a function of", var, "for", type)) +
    theme_minimal()
  
  if (!is.null(fit_new)){
    p <- p + 
      geom_line(data = pred_df, 
        aes(x = .data[[var]], y = hr_new),
        colour = "red") +
      geom_ribbon(data = pred_df,
        aes(x = .data[[var]], ymin = lower_new, ymax = upper_new),
        alpha = 0.2, fill = "red")
  }
  
  if (save == TRUE){
    
    filename = glue("splineplot_{type}_{var}.png")
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
  
  p
}
