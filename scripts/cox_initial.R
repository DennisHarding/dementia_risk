library(survival)
library(survminer)

fit <- coxph(Surv(futime, alldem_bin) ~ 
               age_baseline +
               sex + 
               edu +
               gene_apoe +
               smok_ever # + 
            #   alc_freq +
            #   stroke +
            #   alldm +
             #  hf +
             #  ht +
             #  non_hdl
              ,
             data = df)
ph_test <- cox.zph(fit)

print(ph_test)
ggcoxzph(ph_test, ggtheme = theme_bw())

ggcoxzph(ph_test) +
  theme_bw()
