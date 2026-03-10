library(readr)
library(tidyverse)
library(data.table)


censor_date_new <- as.IDate("2026-01-01")

custom_data <- fread("data/custom_data_1.tsv")

df <- custom_data %>%
  rename(
    ad_hes = ad_hes_df,
    vd_hes = vd_hes_df,
    nonad_hes = nonad_hes_df, #vascular-related dementia
    prs = phylo_score
  ) %>%
  mutate(
    exit_date_ad = pmin(death_date, ad_hes, censor_date_new, na.rm=TRUE),
    exit_date_vd = pmin(death_date, vd_hes, censor_date_new, na.rm=TRUE),
    exit_date_nonad = pmin(death_date, nonad_hes, censor_date_new, na.rm=TRUE),
    exit_date_alldem = pmin(death_date, ad_hes, vd_hes, nonad_hes, censor_date_new, na.rm=TRUE),
    futime_ad = as.numeric(difftime(exit_date_ad, entry_date, unit = "days")) / 365.25,
    futime_vd = as.numeric(difftime(exit_date_vd, entry_date, unit = "days")) / 365.25,
    futime_nonad = as.numeric(difftime(exit_date_nonad, entry_date, unit = "days")) / 365.25,
    futime_alldem = as.numeric(difftime(exit_date_alldem, entry_date, unit = "days")) / 365.25,
    age_ad = age_baseline + as.numeric(difftime(ad_hes, entry_date, unit = "days")) / 365.25,
    age_vd = age_baseline + as.numeric(difftime(vd_hes, entry_date, unit = "days")) / 365.25, 
    age_nonad = age_baseline + as.numeric(difftime(nonad_hes, entry_date, unit = "days")) / 365.25, 
    ad_bin = as.numeric(!is.na(ad_hes)),
    vd_bin = as.numeric(!is.na(vd_hes)),
    nonad_bin = as.numeric(!is.na(nonad_hes)),
    alldem_bin = as.numeric(!is.na(coalesce(age_ad, age_vd, age_nonad))),
    edu = case_when(
      edu_yrs == -3 ~ NA_real_,
      edu_yrs >= 9 ~ 1,
      edu_yrs < 9 ~ 0,
      .default = NA_real_
    ),
    edu_cont = case_when(
      edu_yrs == -3 ~ NA_real_,
      edu_yrs == NA ~ NA_real_,
      .default = edu_yrs
    ),
    smok_ever = case_when(
      smok_ever == "Yes" ~ 1,
      smok_ever == "No" ~ 0,
      .default = NA_real_
    ),
    alc = case_when(
      alc_freq %in% c(1,2,3) ~ 0,
      alc_freq %in% c(4,5,6) ~ 1,
      .default = NA_real_
    ),
    #>find when stroke occured before joining study, is there a better way to do this? allowing for updated variables
    stroke_first = pmin(giga_stroke_hes_df, giga_ih_hes_df, is_hes_df, na.rm=TRUE),
    stroke = case_when(
      stroke_first <= entry_date ~ 1,
      stroke_first > entry_date ~ 0,
      .default = NA_real_
    ),
    ht = case_when(
      ht_hes_df <= entry_date ~ 1,
      ht_hes_df > entry_date ~ 0,
      .default = NA_real_
    ),
    hf = case_when(
      hf_hes_df <= entry_date ~ 1,
      hf_hes_df > entry_date ~ 0,
      .default = NA_real_
    ),
    alldm = case_when(
      alldm_hes_df <= entry_date ~ 1,
      alldm_hes_df > entry_date ~ 0,
      .default = NA_real_
    ),
    non_hdl = chol - hdl,
    # Recommendation from Johan to have e33 as reference
    gene_apoe = relevel(as.factor(gene_apoe), ref = "e33")
  ) %>%
  filter(
    British_white == "Caucasian"
  ) %>%
  as_tibble()

# Create all dfs
df_ad <- df %>%
  mutate(
    exit_date = exit_date_ad,
    futime = futime_ad,
    fail_bin = ad_bin,
    fail_age = age_ad
  ) %>%
  filter(
    futime > 0
  )

df_vd <- df%>%
  mutate(
    exit_date = exit_date_vd,
    futime = futime_vd,
    fail_bin = vd_bin,
    fail_age = age_vd
  ) %>%
  filter(
    futime > 0
  )

df_nonad <- df %>%
  mutate(
    exit_date = exit_date_nonad,
    futime = futime_nonad,
    fail_bin = nonad_bin,
    fail_age = age_nonad
  ) %>%
  filter(
    futime > 0
  )

df_alldem <- df %>%
  mutate(
    exit_date = exit_date_alldem,
    futime = futime_alldem,
    fail_bin = alldem_bin,
  ) %>%
  filter(
    futime > 0
  )
