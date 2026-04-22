fgr <- function(type, data){
  min_sel <- list("AD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "edu"),
                  "VD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "alldm", "is_ih", "ihd", "ht"),
                  "VRD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "edu", "alldm", "ihd", "ht"))
  
  mid_sel <- list("AD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "edu"),
                  "VD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "alldm", "is_ih", "ihd", "ht"),
                  "VRD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "edu", "alldm", "ihd", "ht"))
  
  max_sel <- list("AD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "edu"),
                  "VD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "alldm", "is_ih", "ihd", "ht"),
                  "VRD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "edu", "alldm", "ihd", "ht"))
  
  print(glue("Fitting FGR for {type}"))
  formula <- as.formula(paste0("Hist(futime, fail_cr) ~ ", paste0(cox_selected[[type]], collapse = " + ")))
  fit_fgr <- FGR(formula = formula, data = data, cause = 1)
  print(warnings())
  print(glue("FGR for {type} Fitted"))
  saveRDS(fit_fgr, glue("fits/{type}_fgr.rds"))
  fit_fgr
}
