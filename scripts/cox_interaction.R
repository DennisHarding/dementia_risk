install.packages("survival")
library(survival)
install.packages("survminer")
library(survminer)



fit_main <- coxph(Surv(futime, alldem_bin) ~ 
                    age_baseline + factor(alldm), data = df)


fit <- coxph(Surv(futime, alldem_bin) ~ 
               age_baseline^2 +
               age_baseline +
               factor(sex) + 
               factor(edu) +
               factor(gene_apoe) +
               factor(smok_ever)
             ,
             data = df)
fit2 <- coxph(Surv(futime, alldem_bin) ~ 
               age_baseline +
               factor(sex) + 
               factor(edu) +
               factor(gene_apoe) +
               factor(smok_ever)
             ,
             data = df)



fit_interact <- coxph(Surv(futime, alldem_bin) ~ 
                    age_baseline * factor(alldm), data = df) + age_baseline

anova(fit_main, fit_interact)


AIC(fit, fit2)
