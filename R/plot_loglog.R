plot_loglog <- function(type, data, vars, save = FALSE){
  plots <- map(vars, function(var){
    
    formula <- as.formula(paste("Surv(futime, event = fail_bin) ~", var))
    
    fit <- surv_fit(formula, data = data)
    
    title <- glue("{var} in {type}")
    
    min_fit <- min(fit$time[fit$surv != 1])
    
    max_fit <- max(fit$time)
    
    p <- ggsurvplot(fit = fit,
               fun = "cloglog",
               palette = "Set1",
               legend = "bottom",
               legend.title = title,
               xlab = "Futime",
               ylab = "log(-log(Survival))",
               ggtheme = theme_minimal(),
               xlim = c(min_fit, max_fit),
               censor = FALSE,
               linewidth = 0.5
              )
    
    p$plot
  })
  
  combined <- wrap_plots(plots, ncol = 3)
  
  if (save == TRUE){
    vars_str <- paste0(vars, collapse = "_")
    filename = glue("loglogplot_{type}_{vars_str}.png")
    path = "figures/loglogplot/"
    
    print(glue("Saving {filename} into {path}"))
    
    ggsave(
      filename,
      width = 15, 
      height = 5, 
      plot = combined,
      path = path
    )
  }
  combined
}