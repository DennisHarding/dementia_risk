0, 1, 2, = c("censor", "fail", "competing risk")

cox_selected <- list("AD" = c("sex", "prs", "gene_apoe", "edu", "ht"),
                     "VD" = c("sex", "prs", "gene_apoe", "alldm", "is_ih", "ht"),
                     "VRD" = c("sex", "prs", "gene_apoe", "alldm", "is_ih", "ht"))

# extract 10 percent of dataset in respect to the competing risk variable.
split <- initial_split(df_ad_test, prop = 0.1, strata = fail_cr)

df_train_fgr <- training(split)

hist <- Hist(time = futime, status = fail_cr)


print("Fitting FGR")
formula <- as.formula(paste0("Hist(", age, ", ", bin, ", entry = age_baseline) ~ ", paste0(cox_selected[[pheno]], collapse = " + ")))
fg <- FGR(formula = formula, data=df_train, cause = 1)
print(warnings())
print("FGR Fitted")
saveRDS(fine_gray, paste0("fgr_", tolower(pheno), ".rds"))
print(summary(fine_gray))

save(fg)