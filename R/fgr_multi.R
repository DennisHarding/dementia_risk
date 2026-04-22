fgr_multi <- function(type, data){
  
  selections <- list(
    
    min = list("AD" = c("age_baseline", "gene_apoe", "edu"), # 3 vars
               "VD" = c("age_baseline", "gene_apoe", "alldm", "ihd", "ht"), # 5 vars
               "VRD" = c("age_baseline", "gene_apoe", "edu", "alldm", "ihd", "ht") # 6 vars
    ), 
    
    mid = list("AD" = c("age_baseline", "prs_level", "gene_apoe", "edu", "alldm", "ihs", "ht", "alc"), # 8 vars
               "VD" = c("age_baseline", "sex", "prs_level", "gene_apoe", "edu", "alldm", "is_ih", "ihd", "ht"), # 9 vars
               "VRD" = c("age_baseline", "sex", "prs_level", "gene_apoe", "edu", "alldm", "is_ih", "ihd", "ht", "alc") # 10 vars
    ),
    
    max = list("AD" = c("age_baseline", "sex", "prs_level", "gene_apoe_level", "edu", "smok_ever", "alc", "alldm", "ht", "is_ih", "ihd", "dbp10", "sbp10"), # 13 vars
               "VD" = c("age_baseline", "sex", "prs_level", "gene_apoe_level", "edu", "smok_ever", "alc", "alldm", "ht", "is_ih", "ihd", "dbp10", "sbp10"), # 13 vars
               "VRD" = c("age_baseline", "sex", "prs_level", "gene_apoe_level", "edu", "smok_ever", "alc", "alldm", "ht", "is_ih", "ihd", "dbp10", "sbp10") # 13 vars
    )
    
  )
  
  fit_one <- function(sel_name, vars){
    message(glue("Fitting FGR [{sel_name}] for {type}"))
    
    formula <- as.formula(
      paste0("Hist(futime, fail_cr) ~ ", paste0(vars, collapse = " + "))
    )
    
    fit <- FGR(formula = formula, data = data, cause = 1)
    message(glue("FGR [{sel_name}] for {type} Fitted"))
    path <- glue("fits/{type}_{sel_name}_fgr.rds")
    saveRDS(fit, path)
    fit
  }
  
  fits <- imap(selections, ~ fit_one(.y, .x[[type]]))
  
  fits
}
