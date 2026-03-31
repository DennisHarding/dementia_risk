

lbs_fun <- function(lra, ...) {
  
  fit <- lra$glmnet.fit
  
  L=which(fit$lambda==lra$lambda.min)
  
  ystart <- sort(fit$beta[abs(fit$beta[,L])>0,L])
  labs <- names(ystart)
  r <- range(fit$beta[,100]) # max gap between biggest and smallest coefs at smallest lambda i.e., 100th lambda
  yfin <- seq(r[1],r[2],length=length(ystart))
  
  xstart<- log(lra$lambda.min)
  xfin <- xstart+1
  
  
  text(xfin+0.3,yfin,labels=labs,...)
  segments(xstart,ystart,xfin,yfin)
  
  
}

plot(lra$glmnet.fit,label=F, xvar="lambda", xlim=c(-5.2,0), lwd=2) #xlim, lwd is optional


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

data <- df_ad
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

head(X)
y <- Surv(data$futime, event = data$fail_bin)
complete_cases <- complete.cases(X, y)
X_clean <- model.matrix(~ ., data = X[complete_cases, ])[,-1]
y_clean <- y[complete_cases]
head(X_clean)
head(y_clean)
fit <- glmnet(X_clean, y_clean, family = "cox")

plot(fit, label = TRUE)

cvfit <- cv.glmnet(X_clean, y_clean, family = "cox")

plot(cvfit)

-log(cvfit$lambda.min)
-log(cvfit$lambda.1se)

coef(cvfit, s = "lambda.1se")
coef(cvfit, s = "lambda.min")


