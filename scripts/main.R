#> -----------------------------
#> LOAD LIBRARIES
#> -----------------------------
{
library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(glue)
library(gtable)
library(patchwork)
library(broom)
library(glmnet)
library(forestplot)
library(rsample)
library(riskRegression)
library(prodlim)
library(splines)
library(furrr)
library(pec)
}
#> -----------------------------
#> DATA FORMATTING
#> -----------------------------
{
set.seed(10)

UKB_data <- fread("data/custom_data_1.tsv")

source("R/formatting.R")

df <- format_df(UKB_data)

ad_split <- format_df_ad(df)
vd_split <- format_df_vd(df)
vrd_split <- format_df_vrd(df)

df_ad_train <- training(ad_split)
df_vd_train <- training(vd_split)
df_vrd_train <- training(vrd_split)

df_ad_test <- testing(ad_split)
df_vd_test <- testing(vd_split)
df_vrd_test <- testing(vrd_split)
}
#> -----------------------------
#> Prepare analysis table
#> -----------------------------
{
types_list <- list(
  AD = "AD",
  VD = "VD",
  VRD = "VRD"
)

data_train_list <- list(
  AD_train = df_ad_train,
  VD_train = df_vd_train,
  VRD_train = df_vrd_train
  
)

data_test_list <- list(
  AD_test = df_ad_test,
  VD_test = df_vd_test,
  VRD_test = df_vrd_test
)

analysis_tbl <- tibble(
  type = types_list,
  data_train = data_train_list,
  data_test = data_test_list
)
}
#> -----------------------------
#> Forest plots
#> -----------------------------
{
source("R/plot_forest.R")

analysis_tbl = analysis_tbl %>% 
  mutate(forest_plot = pmap(list(type, data_train), ~ 
                              plot_forest(..1, ..2)))
}
levels(df$prs_fac)
#> -----------------------------
#> DEFINE & FIT COX MODELS
#> -----------------------------

{analysis_tbl <- analysis_tbl %>%
  mutate(
    data_train = map(data_train, ~ filter(.x, !is.na(sbp),
                              !is.na(dbp),
                              !is.na(ht))
    )
  )
analysis_tbl <- analysis_tbl %>% 
  mutate(
    cox_fit = map(data_train, ~ coxph(
      Surv(futime, fail_bin) ~ 
        prs +
        sex + 
        gene_apoe +
        smok_ever +
        alc +
        edu_cont +
        age_baseline +
        dbp10 +
        sbp10 +
        ht +
        ihd +
        is_ih
      ,
      data = .x
    ))
  )

analysis_tbl %>% 
  mutate(cox_summary = map(cox_fit, summary)) %>% 
  pull(cox_summary) %>% 
  walk(print)
# Optional: define candidate cox model:

analysis_tbl <- analysis_tbl %>% 
  mutate(
    cox_fit_new = map(data_train, ~ coxph(
      Surv(futime, fail_bin) ~ 
        I(age_baseline^2) +
        sex + 
        prs +
        gene_apoe +
        edu_cont +
        smok_ever +
        alc +
        ht

      ,
      data_train = .x
    ))
  )


analysis_tbl %>% 
  mutate(cox_summary_new = map(cox_fit_new, summary)) %>% 
  pull(cox_summary_new) %>% 
  walk(print)
}
#> -----------------------------
#> fgr run 
#> -----------------------------

cox_selected <- list("AD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "edu"),
                     "VD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "alldm", "is_ih", "ihd", "ht"),
                     "VRD" = c("age_gr", "sex", "prs_fac", "gene_apoe", "edu", "alldm", "ihd", "ht"))

formula <- as.formula(paste0("Hist(futime, fail_cr) ~ ", paste0(cox_selected[["AD"]], collapse = " + ")))
fit_fgr <- FGR(formula = formula, data = df_ad_train, cause = 1)
saveRDS(fit_fgr, glue("fits/{type}_fgr.rds"))
fit_fgr
df %>% 
  select(cox_selected[["AD"]]) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  cor(use = "complete.obs")

{
source("R/fgr.R")
all_deps <- tools::package_dependencies("riskRegression", recursive = TRUE)$riskRegression

plan(multisession, workers = 3)

analysis_tbl <- analysis_tbl %>%
  mutate(fgr_fit = future_pmap(
    list(type, data_train),
    ~ fgr(..1, ..2),
    .options = furrr_options(
      packages = c("Matrix", "riskRegression", all_deps, "rsample", "glue", "dplyr"),
      seed = TRUE
    )
  ))
}
#> -----------------------------
#> Risk Charts
#> -----------------------------

AD.fit <- readRDS("fits/AD_fgr.rds")
summary(AD.fit)

VD.fit <- readRDS("fits/VD_fgr.rds")
summary(VD.fit)

VRD.fit <- readRDS("fits/VRD_fgr.rds")
summary(VRD.fit)


## 10 year Risk estimates

# r setup and compute 10 year risk estimates
# Setup array for results

# loop - sex
# loop - age_group
# loop - APOE_genotype
# loop - prs_fac
# loop - Diabetes status
# loop - Smoking_n_e
# loop - Education_gr

sex <- as.factor(c(0,1))
age_gr <- as.factor(c("40-50", "50-60", "60-70","70-80"))
gene_apoe <- relevel(as.factor(c("e33","e22","e32","e42","e43","e44")),ref="e33")
prs_fac <- relevel(as.factor(c("1","2","3","4","5")),ref = "3")
alldm <- as.factor(c(0,1))
ihd <- as.factor(c(0,1))
is_ih <- as.factor(c(0,1))
ht <- as.factor(c(0,1))
smoker <- as.factor(c(0,1))
edu <- as.factor(c(0,1))

m1.AD <- CJ(age_gr, sex, prs_fac, gene_apoe, edu)
m1.VD <- CJ(age_gr, sex, prs_fac, gene_apoe, alldm, is_ih, ihd, ht)
m1.VRD <- CJ(age_gr, sex, prs_fac, gene_apoe, edu, alldm, ihd, ht)


rm(sex,age_gr,gene_apoe,prs_fac,alldm,ihd,is_ih,ht,smoker,edu)

# Risk prediction times 5,10,15,20 years

pred_ad <- as.data.frame(predictEventProb(AD.fit, times=c(5,10,15,20), newdata=m1.AD, cause=1))
colnames(pred_ad) <- c("ad_pred_05","ad_pred_10","ad_pred_15","ad_pred_20")

pred_vd <- as.data.frame(predictEventProb(VD.fit, times=c(5,10,15,20), newdata=m1.VD, cause=1))
colnames(pred_vd) <- c("vd_pred_05","vd_pred_10","vd_pred_15","vd_pred_20")

pred_vrd <- as.data.frame(predictEventProb(VRD.fit, times=c(5,10,15,20), newdata=m1.VRD, cause=1))
colnames(pred_vrd) <- c("vrd_pred_05","vrd_pred_10","vrd_pred_15","vrd_pred_20")



ad_df <- cbind(m1.AD,pred_ad)
vd_df <- cbind(m1.VD,pred_vd)
vrd_df <- cbind(m1.VRD,pred_vrd)

rm(pred_ad,pred_vd,pred_vrd,m1.AD,m1.VD,m1.VRD)

# Set up variables for risk plots

age_gr.labs <- c("50-59", "60-69", "70-79",">80")
names(age_gr.labs) <- c("1", "2", "3","4")

smoker.labs <- c("No smoking", "Smoking")
names(smoker.labs) <- c("0", "1")

edu.labs <- c("Education >9y", "Education <9y")
names(edu.labs) <- c("1", "0")

sex.labs <- c("Male", "Female")
names(sex.labs) <- c("1, 0")

dm.labs <- c("No diabetes", "Diabetes")
names(dm.labs) <- c("0", "1")

cols <- c("0"="#00994C","1"="#66CC00","2"="#FFFF33" ,"3"="#FF9933", "4"="#FF0000","5"="#CC0000","6"="#994C00")

# adheimer's disease risk charts

ad_df <- ad_df %>% mutate(apoe_s = case_when(
  gene_apoe == "e22" ~ 1,
  gene_apoe == "e32" ~ 2,
  gene_apoe == "e33" ~ 3,
  gene_apoe == "e42" ~ 4,
  gene_apoe == "e43" ~ 5,
  gene_apoe == "e44" ~ 6)
)

ad_df$apoe_s<- as.factor(ad_df$apoe_s)

ad_df <- ad_df %>% mutate(col = case_when(
  ad_pred_10 < 0.01 ~ "0",
  ad_pred_10 >= 0.01 & ad_pred_10 < 0.05 ~ "1",
  ad_pred_10 >= 0.05 & ad_pred_10 < 0.10 ~ "2",
  ad_pred_10 >= 0.10 & ad_pred_10 < 0.15 ~ "3",
  ad_pred_10 >= 0.15 & ad_pred_10 < 0.20 ~ "4",
  ad_pred_10 >= 0.20 & ad_pred_10 < 0.30 ~ "5",
  ad_pred_10 >= 0.30 ~ "6")
)

# order age_gr for plot 
summary(ad_df$age_gr)
ad_df$age_gr <- factor(ad_df$age_gr, levels = c("40-50", "50-60", "60-70","70-80"))
ad_df$prs_fac<- factor(ad_df$prs_fac, levels = c("1","2","3","4","5"))
ad_df$prs_fac

ad_df

pd_ad <- ad_df[age_gr != 0,]
titlename <- paste0("dharding - alzheimer 10 year absolute risk - sex=" ,i, " edu=" ,j)
pd_ad
ggplot(pd_ad,aes(x = apoe_s, y = prs_fac, fill = col))  +
  geom_tile() +
  geom_text(aes(label = lapply(ad_pred_10*100.0, as.integer)),fontface = "bold") +
  scale_fill_manual(values=cols)+
  labs(x="APOE genotype", y="Allele score") + 
  scale_x_discrete( labels=c(expression(paste(epsilon,"22")),
                             expression(paste(epsilon,"32")),
                             expression(paste(epsilon,"33")),
                             expression(paste(epsilon,"42")),
                             expression(paste(epsilon,"43")),
                             expression(paste(epsilon,"44"))
                             )) + 
  scale_y_discrete(labels=c(
    "1" = "0-20%",
    "2" = "20-40%",
    "3" = "40-60%",
    "4" = "60-80%",
    "5" = "80-100%"
  )) + 
  theme(axis.ticks = element_blank(), 
        panel.background = element_blank(),
        legend.position = "none",
        axis.title=element_text(size=12,face="bold")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0 , b = 0, l = 0))) +
  theme(axis.text=element_text(size=8,face="italic"), strip.background = element_blank()) +
  facet_grid(age_gr ~ sex + edu , labeller = labeller(age_gr = age_gr.labs, sex = sex.labs, edu = edu.labs)) +
  theme(strip.text.x = element_text( size = 10, face = "bold"),
        strip.text.y = element_text( size = 10, face = "bold")) +
  ggtitle(titlename) +
  coord_fixed()




for (i in 0:1) {
  for (j in 0:1) {
    
    pd_ad <- ad_df[sex==i & edu==j & age_gr != 0,]
    
    titlename <- paste0("dharding - alzheimer 10 year absolute risk - sex=" ,i, " edu=" ,j)
    
    ggplot(pd_ad,aes(x = apoe_s, y = prs_fac, fill = col))  +
      geom_tile() +
      geom_text(aes(label = lapply(ad_pred_10*100.0, as.integer)),fontface = "bold") +
      scale_fill_manual(values=cols)+
      labs(x="APOE genotype", y="Allele score") + 
      scale_x_discrete( labels=c(expression(paste(epsilon,"22")),
                                 expression(paste(epsilon,"32")),
                                 expression(paste(epsilon,"33")),
                                 expression(paste(epsilon,"42")),
                                 expression(paste(epsilon,"43")),
                                 expression(paste(epsilon,"44"))
                                 )) + 
      scale_y_discrete( labels=c("0-20%","20-40%","40-60%","60-80%", "80-100%")) + 
      theme(axis.ticks = element_blank(), panel.background = element_blank(),legend.position = "none",      axis.title=element_text(size=12,face="bold")) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
      theme(axis.title.x = element_text(margin = margin(t = 10, r = 0 , b = 0, l = 0))) +
      theme(axis.text=element_text(size=8,face="italic"), strip.background = element_blank()) +
      facet_grid(age_gr ~ sex + edu , labeller = labeller(age_gr = age_gr.labs, sex = sex.labs, edu = edu.labs)) +
      theme(strip.text.x = element_text( size = 10, face = "bold"),
            strip.text.y = element_text( size = 10, face = "bold")) +
      ggtitle(titlename) +
      coord_fixed()
    
    
    plotname <- paste0("./Results/dharding_EADB2022_Riskplot_10y_ad_sex",i,"_edu",j,".pdf")
    
    ggsave(plotname,device="pdf",width=20,dpi=300,units="cm")
    
  }
}

# Vascular dementia risk charts

vd_df <- vd_df %>% mutate(apoe_s = case_when(
  gene_apoe == "e22" ~ 1,
  gene_apoe == "e32" ~ 2,
  gene_apoe == "e33" ~ 3,
  gene_apoe == "e42" ~ 4,
  gene_apoe == "e43" ~ 5,
  gene_apoe == "e44" ~ 6)
)

vd_df$apoe_s<- as.factor(vd_df$apoe_s)

vd_df <- vd_df %>% mutate(col = case_when(
  vd_pred_10 < 0.01 ~ "0",
  vd_pred_10 >= 0.01 & vd_pred_10 < 0.05 ~ "1",
  vd_pred_10 >= 0.05 & vd_pred_10 < 0.10 ~ "2",
  vd_pred_10 >= 0.10 & vd_pred_10 < 0.15 ~ "3",
  vd_pred_10 >= 0.15 & vd_pred_10 < 0.20 ~ "4",
  vd_pred_10 >= 0.20 & vd_pred_10 < 0.30 ~ "5",
  vd_pred_10 >= 0.30 ~ "6")
)

# order age_gr for plot 
summary(vd_df$age_gr)
vd_df$age_gr <- factor(vd_df$age_gr, levels = c(4, 3, 2, 1,0))


for (i in 0:1) {
  for (j in 0:1) {
    
    pd_vd <- vd_df[sex==i & edu==j & age_gr != 0,]
    
    titlename <- paste0("dharding - vascular dementia 10 year absolute risk - sex=" ,i, " edu=" ,j)
    
    ggplot(pd_vd,aes(x = apoe_s, y = was2022_gr, fill = col))  +
      geom_tile() +
      geom_text(aes(label = lapply(vd_pred_10*100.0, as.integer)),fontface = "bold") +
      scale_fill_manual(values=cols)+
      labs(x="APOE genotype", y="Allele score") + 
      scale_x_discrete( labels=c(expression(paste(epsilon,"22")),expression(paste(epsilon,"32")),expression(paste(epsilon,"33")),expression(paste(epsilon,"42")),expression(paste(epsilon,"43")),expression(paste(epsilon,"44")))) + 
      scale_y_discrete( labels=c("0-25%","25-50%","50-75%","75-100%")) + 
      theme(axis.ticks = element_blank(), panel.background = element_blank(),legend.position = "none",      axis.title=element_text(size=12,face="bold")) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
      theme(axis.title.x = element_text(margin = margin(t = 10, r = 0 , b = 0, l = 0))) +
      theme(axis.text=element_text(size=8,face="italic"), strip.background = element_blank()) +
      facet_grid(age_gr ~ dm + smoker , labeller = labeller(age_gr = age_gr.labs, smoker = smoker.labs, dm = dm.labs)) +
      theme(strip.text.x = element_text( size = 10, face = "bold"),
            strip.text.y = element_text( size = 10, face = "bold")) +
      ggtitle(titlename) +
      coord_fixed()
    
    
    plotname <- paste0("./Results/dharding_EADB2022_Riskplot_10y_vd_sex",i,"_edu",j,".pdf")
    
    ggsave(plotname,device="pdf",width=20,dpi=300,units="cm")
    
  }
}

# vdcular related disease risk charts


vrd_df <- vrd_df %>% mutate(apoe_s = case_when(
  gene_apoe == "e22" ~ 1,
  gene_apoe == "e32" ~ 2,
  gene_apoe == "e33" ~ 3,
  gene_apoe == "e42" ~ 4,
  gene_apoe == "e43" ~ 5,
  gene_apoe == "e44" ~ 6)
)

vrd_df$apoe_s<- as.factor(vrd_df$apoe_s)

vrd_df <- vrd_df %>% mutate(col = case_when(
  vrd_pred_10 < 0.01 ~ "0",
  vrd_pred_10 >= 0.01 & vrd_pred_10 < 0.05 ~ "1",
  vrd_pred_10 >= 0.05 & vrd_pred_10 < 0.10 ~ "2",
  vrd_pred_10 >= 0.10 & vrd_pred_10 < 0.15 ~ "3",
  vrd_pred_10 >= 0.15 & vrd_pred_10 < 0.20 ~ "4",
  vrd_pred_10 >= 0.20 & vrd_pred_10 < 0.30 ~ "5",
  vrd_pred_10 >= 0.30 ~ "6")
)

# order age_gr for plot 
summary(vrd_df$age_gr)
vrd_df$age_gr <- factor(vrd_df$age_gr, levels = c(4, 3, 2, 1,0))


for (i in 0:1) {
  for (j in 0:1) {
    
    pd_vrd <- vrd_df[sex==i & edu==j & age_gr != 0,]
    
    titlename <- paste0("dharding - vascular related dementia 10 year absolute risk - sex=" ,i, " edu=" ,j)
    
    ggplot(pd_vrd,aes(x = apoe_s, y = was2022_gr, fill = col))  +
      geom_tile() +
      geom_text(aes(label = lapply(vrd_pred_10*100.0, as.integer)),fontface = "bold") +
      scale_fill_manual(values=cols)+
      labs(x="APOE genotype", y="Allele score") + 
      scale_x_discrete( labels=c(expression(paste(epsilon,"22")),expression(paste(epsilon,"32")),expression(paste(epsilon,"33")),expression(paste(epsilon,"42")),expression(paste(epsilon,"43")),expression(paste(epsilon,"44")))) + 
      scale_y_discrete( labels=c("0-25%","25-50%","50-75%","75-100%")) + 
      theme(axis.ticks = element_blank(), panel.background = element_blank(),legend.position = "none",      axis.title=element_text(size=12,face="bold")) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
      theme(axis.title.x = element_text(margin = margin(t = 10, r = 0 , b = 0, l = 0))) +
      theme(axis.text=element_text(size=8,face="italic"), strip.background = element_blank()) +
      facet_grid(age_gr ~ dm + smoker , labeller = labeller(age_gr = age_gr.labs, smoker = smoker.labs, dm = dm.labs)) +
      theme(strip.text.x = element_text( size = 10, face = "bold"),
            strip.text.y = element_text( size = 10, face = "bold")) +
      ggtitle(titlename) +
      coord_fixed()
    
    
    plotname <- paste0("./Results/dharding_EADB2022_Riskplot_10y_vrd_sex",i,"_edu",j,".pdf")
    
    ggsave(plotname,device="pdf",width=20,dpi=300,units="cm")
    
  }
}


## Generate 10 year risk-ratios
#```{r generate risk ratios}
# Women, no diabetes, never smoking, e33, 3rd allele_gr, age_gr 70-80, Normal education is selected as reference

# AD
ref_m <- filter(ad_df, sex==0 & dm==0 & smoker == 0 &  edu == 0 & gene_apoe =="e33" & was2022_gr == 2 & age_gr == 3)
ref <- as.double(ref_m$ad_pred_10)


ad_df <-ad_df %>% mutate( ad_rr_10 = ad_pred_10 / ref)

# VD
ref_m <- filter(vd_df, sex==0 & dm==0 & smoker == 0 &  edu == 0 & gene_apoe =="e33" & was2022_gr == 2 & age_gr == 3)
ref <- as.double(ref_m$vd_pred_10)


vd_df <-vd_df %>% mutate( vd_rr_10 = vd_pred_10 / ref)

# VRD 
ref_m <- filter(vrd_df, sex==0 & dm==0 & smoker == 0 &  edu == 0 & gene_apoe =="e33" & was2022_gr == 2 & age_gr == 3)
ref <- as.double(ref_m$vrd_pred_10)


vrd_df <- vrd_df %>% mutate( vrd_rr_10 = vrd_pred_10 / ref)


#> -----------------------------
#> Anova test
#> -----------------------------
{
analysis_tbl <- analysis_tbl %>%
  mutate(
    anova_test = pmap(list(cox_fit, cox_fit_new), ~ 
      anova(..1,..2)
    )
  )

# Print anova tests
analysis_tbl %>% pull(anova_test) %>% walk(print)
}
#> -----------------------------
#> AIC test
#> -----------------------------
{
analysis_tbl <- analysis_tbl %>%
  mutate(
    AIC_test = pmap(list(cox_fit, cox_fit_new), ~ 
      AIC(..1,..2)
    )
  )

# Print AIC tests
analysis_tbl %>% pull(AIC_test) %>% walk(print)

}
#> -----------------------------
#> SPLINES
#> -----------------------------
{
source("R/plot_spline_dfs.R")
spline_var <- "age_baseline"
  
spline_vars <- c(    
  "prs", # continuous
  "sex", 
  "gene_apoe",
  "smok_ever",
  "alc",
  "edu_cont", # continuous
  "age_baseline", # continuous
  "dbp10", # continuous
  "sbp10", # continuous
  "ht",
  "is_ih",
  "ihd"
)

# IF VARS IS MODIFIED: EDIT create_pred_df() IN R/plot_spline_dfs.R

analysis_tbl <- analysis_tbl %>%
  mutate(
    multispline_plot = pmap(
      list(type, data_train), ~
        plot_spline_dfs(..1, ..2, spline_vars, spline_var, 1, 3, save = TRUE)
        )
  )

}

{
source("R/plot_spline.R")

analysis_tbl <- analysis_tbl %>% 
  mutate(
    cox_fit_spline = map(data_train, ~ coxph(
      Surv(futime, fail_bin) ~ 
        age_baseline +
        sex + 
        prs +
        gene_apoe +
        edu_cont +
        smok_ever +
        alc +
        dbp10 +
        sbp10
      ,
      data = .x
    ))
  )
# Generate reference prediction data_train
analysis_tbl <- analysis_tbl %>% 
  mutate(
    pred_df = map(data_train, ~ tibble(
      prs = mean(.x$prs, na.rm = TRUE),
      sex = as.factor(0),
      edu = as.factor(0),
      gene_apoe = factor("e33", levels = levels(.x$gene_apoe)),
      smok_ever = as.factor(0),
      alc = as.factor(0),
      edu_cont = mean(.x$edu_cont, na.rm = TRUE),
      stroke = as.factor(0),
      age_baseline = mean(.x$age_baseline, na.rm = TRUE),
      dbp10 = as.factor(0),
      sbp10 = seq(min(.x$sbp, na.rm = TRUE),
                max(.x$sbp, na.rm = TRUE),
                length.out = 200)
    ))
  )

# Produce spline plots
analysis_tbl <- analysis_tbl %>%
  mutate(
    spline_plot = pmap(list(type, cox_fit_spline, data_train, pred_df),
               ~ plot_spline(..1, ..2, ..3, ..4, 
                             var = "sbp10",
                             save = TRUE))
  )

# PRINT SPLINE PLOTS
analysis_tbl %>% pull(spline_plot) %>% walk(print)

# # Produce spline plots WITH COMPARISON
# analysis_tbl <- analysis_tbl %>%
#   mutate(
#     spline_plot_vs = pmap(list(type, cox_fit, data_train, pred_df, cox_fit_new),
#                ~ plot_spline(..1, ..2, ..3, ..4, ..5, 
#                              var_name = "age_baseline",
#                              save = FALSE)
#     )
#   )
# 
# 
# # Produce spline plots
# analysis_tbl %>% pull(spline_plot_vs) %>% walk(print)
}
#> -----------------------------
#> Proportional Hazards tests - Schoenfeld residuals
#> -----------------------------
{source("R/plot_sfr.R")

# Tests
analysis_tbl <- analysis_tbl %>%
  mutate(
    ph_test = map(cox_fit, ~
      cox.zph(.x)
    )
  )


# Plots
analysis_tbl <- analysis_tbl %>%
  mutate(
    ph_plot = pmap(list(type, ph_test), ~ 
      plot_sfr(..1, ..2, save = FALSE))
  )

# PRINT PH PLOTS
analysis_tbl %>% pull(ph_plot) %>% walk(print)
}
#> -----------------------------
#> Proportional Hazards tests - Log-Log plots
#> -----------------------------
{source("R/plot_loglog.R")

analysis_tbl <- analysis_tbl %>%
  mutate(surv_curvs = pmap(list(type, data_train), ~ {
    
    # SPECIFY vars of interest -->
    vars = c("ht", "ihd", "sex")
    
    plot_loglog(..1, ..2, vars, save = TRUE)
  }))

analysis_tbl %>% pull(surv_curvs) %>% walk(print)
}
#> -----------------------------
#> Feature selection - lasso cox
#> -----------------------------
{
source("R/plot_lasso.R")

lasso_vars <- c(
  "age_baseline", 
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
  "sbp10")

lasso_selected <- list("AD" = c("age_baseline", "sex"),
                     "VD" = c("age_baseline", "sex"),
                     "VRD" = c("age_baseline", "sex"))

analysis_tbl <- analysis_tbl %>%
  mutate(lasoo_plot = pmap(list(type, data_train), ~ {
    plot_lasso(..1, ..2, lasso_vars, lasso_selected)
  }))

}
#> -----------------------------
#> Generate Summary Table
#> -----------------------------
{
library(gtsummary)

#remember to select variables to include in summary table

summary_table <- 
  df %>%
  select(
    "age_baseline",
    "sex",
    "edu_cont",
    "edu",
    "gene_apoe",
    "alldem_bin",
    "ad_bin",
    "vd_bin",
    "vrd_bin",
    "prs_fac",
    "alldm",
    "sbp",
    "dbp"
  ) %>%
  mutate(
    sex = factor(sex,
                 levels = c(0, 1),
                 labels = c("Female", "Male"))
  ) %>%
  tbl_summary(
    by = sex,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  )

summary_table
}
