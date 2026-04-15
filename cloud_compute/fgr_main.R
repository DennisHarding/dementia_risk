#> -----------------------------
#> INSTALL AND LOAD LIBRARIES
#> -----------------------------

.libPaths(c("R_libs", .libPaths()))

install.packages(c(
  "data.table",
  "dplyr",
  "tibble", 
  "purrr",
  "rsample", 
  "glue", 
  "prodlim",
  "riskRegression"
))

library(data.table)
library(dplyr)
library(tibble)
library(purrr)
library(rsample)
library(glue)
library(prodlim)
library(riskRegression)

#> -----------------------------
#> DATA FORMATTING
#> -----------------------------

set.seed(10)

UKB_data <- fread("custom_data_1.tsv")

source("formatting.R")

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

#> -----------------------------
#> Prepare analysis table
#> -----------------------------

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

#> -----------------------------
#> Fine-Gray test
#> -----------------------------

source("fgr.R")

analysis_tbl <- analysis_tbl %>%
  mutate(fgr_fit = pmap(list(type, data_train), ~
                          fgr(..1, ..2)
  ))



