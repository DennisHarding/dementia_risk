# A function that generates a forest plot based on a coxph object

plot_forest <- function(type, fit){
  # Create a table
  hr_table <- as.data.frame(exp(cbind(HR = coef(fit),
                               Lower = confint(fit)[, 1],
                               Upper = confint(fit)[, 2])))
  # Assign rownames based on the creation of the table
  hr_table$variable <- rownames(hr_table)
  
  # Add pvalues
  hr_table$pvalue <- summary(fit)$coefficients[, "Pr(>|z|)"]
   
  hr_table$p_text <- sprintf("%.2f (%.2f–%.2f)", 
                             hr_table$HR, 
                             hr_table$Lower, 
                             hr_table$Upper)
  hr_table$se <- (log(hr_table$Upper) - log(hr_table$Lower))/(2*1.96)
  
  names = names(coef(fit))
  title = glue("Cox regression results for: {type}")
  
  forestplot(df = hr_table,
             name = names,
             estimate = HR,
             se = se,
             pvalue = pvalue,
             psignif = 0.05,
             xlab= "Hazard Ratio (95% CI)",
             title = title
             ) + geom_vline(xintercept = 1, linetype = "solid", colour = "red")
  
}
