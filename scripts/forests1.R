library(survival)
remotes::install_github("NightingaleHealth/ggforestplot")
library(ggforestplot)

#fit a model to chosen covariates

fit <- coxph(Surv(futime, alldem_bin) ~ 
  age_baseline + 
  factor(sex) + 
  factor(edu) + 
  factor(gene_apoe) +
  factor(smok_ever), 
  data = df)
summary(fit)

# Create a table
hr_table <- as.data.frame(exp(cbind(HR = coef(fit),
                             Lower = confint(fit)[, 1],
                             Upper = confint(fit)[, 2])))
# Assign rownames based on the creation of the table
hr_table$variable <- rownames(hr_table)

# Add pvalues
hr_table$pvalue <- summary(fit)$coefficients[, "Pr(>|z|)"]


hr_table$label <- c("Age_baseline", 
                    "sex", 
                    "edu", 
                    "gene_apoe e23", 
                    "gene_apoe e24", 
                    "gene_apoe e33", 
                    "gene_apoe e34", 
                    "gene_apoe e44", 
                    "smok_ever")
 
hr_table$p_text <- sprintf("%.2f (%.2f–%.2f)", 
                           hr_table$HR, 
                           hr_table$Lower, 
                           hr_table$Upper)
hr_table$se <- (log(hr_table$Upper) - log(hr_table$Lower))/(2*1.96)

glimpse(hr_table)

p <- forestplot(df = hr_table,
           name = label,
           estimate = HR,
           se = se,
           pvalue = pvalue,
           psignif = 0.05,
           ci
           xlab= "Hazard Ratio (95% CI)",
           title = "Cox Regression Results",
           ref_line = 1
           ) + geom_vline(xintercept = 1, linetype = "dashed")

