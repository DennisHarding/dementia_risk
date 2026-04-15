legends <- list(
  "sex" = "Sex", 
  "alldm" = "Diabetes", 
  "smok_ever" = "Smoking", 
  "edu" = "Education", 
  "alc" = "Alcohol intake", 
  "is_ih" = "is_ih Stroke",
  "ht" = "Hypertension", 
  "dbp10" = "Diastolic BP 10",
  "sbp10" = "Systolic BP 10",
  "ihd" = "Ischemic Heart Disease",
  "prs_fac" = "Polygenetic Risk Score"
)

data <- df_ad_train
type = "AD"

X <- data %>%
  select(c("age_baseline", 
           "prs_level", 
           "sex", 
           "gene_apoe_level", 
           "alldm", 
           "smok_ever", 
           "ht", 
           "edu", 
           "alc", 
           "is_ih", 
           "ihd", 
           "dbp10", 
           "sbp10"))

y <- Surv(data$futime, event = data$fail_bin)
complete_cases <- complete.cases(X, y)
X_clean <- model.matrix(~ ., data = X[complete_cases, ])[,-1]
y_clean <- y[complete_cases]


cvfit <- cv.glmnet(X_clean, y_clean, family = "cox", nfolds = 10)

plot(cvfit)

-log(cvfit$lambda.min)
-log(cvfit$lambda.1se)

coef(cvfit, s = "lambda.1se")
coef(cvfit, s = "lambda.min")


