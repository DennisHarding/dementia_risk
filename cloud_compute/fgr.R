fgr <- function(type, data){
  cox_selected <- list("AD" = c("age_baseline", "sex", "prs", "gene_apoe", "edu"),
                       "VD" = c("age_baseline", "sex", "prs", "gene_apoe", "alldm", "is_ih", "ihd", "ht"),
                       "VRD" = c("age_baseline", "sex", "prs", "gene_apoe", "edu", "alldm", "ihd", "ht"))
  
  # extract 10 percent of dataset in respect to the competing risk variable.
  split <- initial_split(data, prop = 0.1, strata = fail_cr)
  
  df_train_fgr <- training(split)
  
  print(glue("Fitting FGR for {type}"))
  formula <- as.formula(paste0("Hist( futime, fail_cr) ~ ", paste0(cox_selected[[type]], collapse = " + ")))
  fit_fgr <- FGR(formula = formula, data = df_train_fgr, cause = 1)
  print(warnings())
  print(glue("FGR for {type} Fitted"))
  saveRDS(fit_fgr, glue("{type}_fgr.rds"))
  print(summary(fit_fgr))
  fit_fgr
}
