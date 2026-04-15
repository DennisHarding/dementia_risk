#> ------------------------------------
#> FUNCTIONS FOR DATA FORMATTING
#> ------------------------------------

format_df <- function(data = data){
  # Update sensor date
  censor_date_new <- as.IDate("2026-01-01")

  df <- data %>%
    # Shorten variable names
    rename(
      age_int = age_baseline,
      ad_hes = ad_hes_df,
      vd_hes = vd_hes_df,
      vrd_hes = nonad_hes_df,   # vascular-related dementia
      prs = phylo_score
    ) %>%
    mutate(
      
      # Define exit dates (earliest of event, death, or censoring)
      exit_date_ad = pmin(death_date, ad_hes, censor_date_new, na.rm=TRUE),
      exit_date_vd = pmin(death_date, vd_hes, censor_date_new, na.rm=TRUE),
      exit_date_vrd = pmin(death_date, vrd_hes, censor_date_new, na.rm=TRUE),
      exit_date_alldem = pmin(death_date, ad_hes, vd_hes, vrd_hes, censor_date_new, na.rm=TRUE),
      
      # Follow-up time (years)
      age_baseline = as.numeric(difftime(entry_date, birth_date, unit = "days")) / 365.25,
      futime_ad = as.numeric(difftime(exit_date_ad, entry_date, unit = "days")) / 365.25,
      futime_vd = as.numeric(difftime(exit_date_vd, entry_date, unit = "days")) / 365.25,
      futime_vrd = as.numeric(difftime(exit_date_vrd, entry_date, unit = "days")) / 365.25,
      futime_alldem = as.numeric(difftime(exit_date_alldem, entry_date, unit = "days")) / 365.25,
      
      # Age at event
      age_ad = age_baseline + as.numeric(difftime(ad_hes, entry_date, unit = "days")) / 365.25,
      age_vd = age_baseline + as.numeric(difftime(vd_hes, entry_date, unit = "days")) / 365.25, 
      age_vrd = age_baseline + as.numeric(difftime(vrd_hes, entry_date, unit = "days")) / 365.25, 
      
      # Event indicators
      ad_bin = as.numeric(!is.na(ad_hes)),
      vd_bin = as.numeric(!is.na(vd_hes)),
      vrd_bin = as.numeric(!is.na(vrd_hes)),
      alldem_bin = as.numeric(!is.na(coalesce(age_ad, age_vd, age_vrd))),
      
      # Sex variable
      sex = factor(sex),
      
      # Education variables
      edu = case_when(
        edu_yrs == -3 ~ NA_real_,   # missing code
        edu_yrs >= 9 ~ 1,
        edu_yrs < 9 ~ 0,
        .default = NA_real_
      ) %>%
        as.factor(),
      edu_cont = case_when(
        edu_yrs == -3 ~ NA_real_,
        edu_yrs == NA ~ NA_real_,
        .default = edu_yrs
      ),
      
      # Smoking (ever vs never)
      smok_ever = case_when(
        smok_ever == "Yes" ~ 1,
        smok_ever == "No" ~ 0,
        .default = NA_real_
      ) %>%
        as.factor(),
      
      # Alcohol frequency grouped
      alc = case_when(
        alc_freq %in% c(1,2,3) ~ 0,
        alc_freq %in% c(4,5,6) ~ 1,
        .default = NA_real_
      ) %>%
        as.factor(),
      
      # Earliest stroke record
      stroke_first = pmin(giga_stroke_hes_df, 
                          giga_ih_hes_df, 
                          is_hes_df, na.rm=TRUE
      ),
      
      # Baseline disease indicators (present before entry)
      allstroke = case_when(
        stroke_first <= entry_date ~ 1,
        stroke_first > entry_date ~ 0,
        .default = 0
      ) %>%
        as.factor(),
      is = case_when(
        is_hes_df <= entry_date ~ 1,
        is_hes_df > entry_date ~ 0,
        .default = 0
      ) %>%
        as.factor(),
      ih = case_when(
        giga_ih_hes_df <= entry_date ~ 1,
        giga_ih_hes_df > entry_date ~ 0,
        .default = 0
      ) %>%
        as.factor(),
      is_ih = case_when(
        is_hes_df <= entry_date | giga_ih_hes_df <= entry_date ~ 1,
        .default = 0
      ) %>%
        as.factor(),
      ihd = case_when(
        ihd_hes_df <= entry_date ~ 1,
        ihd_hes_df > entry_date ~ 0,
        .default = 0
      ) %>%
        as.factor(),
      ht = case_when(
        ht_hes_df <= entry_date ~ 1,
        ht_hes_df > entry_date ~ 0,
        .default = 0
      ) %>%
        as.factor(),
      sbp10 = sbp / 10,
      dbp10 = dbp / 10,
      bpdiff10 = (sbp - dbp) / 10,
      bpadd10 = (sbp + dbp) / 10,
      hf = case_when(
        hf_hes_df <= entry_date ~ 1,
        hf_hes_df > entry_date ~ 0,
        .default = 0
      ) %>%
        as.factor(),
      alldm = case_when(
        alldm_hes_df <= entry_date ~ 1,
        alldm_hes_df > entry_date ~ 0,
        .default = 0
      ) %>%
        as.factor(),
      
      # Derived lipid measure
      non_hdl = chol - hdl,
      non_hdl_bin = case_when(
        non_hdl > median(non_hdl) ~ 1,
        non_hdl <= median(non_hdl) ~ 0,
        .default = NA_real_
      ) %>%
        as.factor(),
      
      prs_fac = factor(ntile(prs, 5), levels = 1:5),
      prs_fac = relevel(prs_fac, ref = "3"),
      prs_level = as.numeric(prs_fac),
      # Set APOE reference genotype
      gene_apoe = relevel(as.factor(gene_apoe), ref = "e33"),
      gene_apoe_level = case_when(
        gene_apoe == "e22" ~ 1,
        gene_apoe == "e23" ~ 2,
        gene_apoe == "e33" ~ 3,
        gene_apoe == "e24" ~ 4,
        gene_apoe == "e34" ~ 5,
        gene_apoe == "e44" ~ 6,
        .default = NA_real_
      ),
      death = !is.na(death_date)
      ) %>%
    
    # Restrict analysis population
    filter(
      British_white == "Caucasian"
    ) %>%
    
    # Ensure tibble output
    as_tibble()
}

# Functions for formatting different data frames

# Alzheimer's Disease
format_df_ad <- function(df = df){ 
  df <- df %>%
    mutate(
      exit_date = exit_date_ad,
      futime = futime_ad,
      fail_bin = ad_bin,
      fail_age = age_ad,
      fail_cr = ifelse(death == TRUE & fail_bin == 0, 2, fail_bin)
    ) %>%
    filter(
      futime > 0
    )
  
  rsample::initial_split(df, prop = 0.8, strata = fail_cr)
}

# Vascular Dementia
format_df_vd <- function(df = df){
  df <- df %>%
    mutate(
      exit_date = exit_date_vd,
      futime = futime_vd,
      fail_bin = vd_bin,
      fail_age = age_vd,
      fail_cr = ifelse(death == TRUE & fail_bin == 0, 2, fail_bin)
    ) %>%
    filter(
      futime > 0
    )
  rsample::initial_split(df, prop = 0.8, strata = fail_cr)
}

# Vascular-Related Dementia
format_df_vrd <- function(df = df){
  df <- df %>%
    mutate(
      exit_date = exit_date_vrd,
      futime = futime_vrd,
      fail_bin = vrd_bin,
      fail_age = age_vrd,
      fail_cr = ifelse(death == TRUE & fail_bin == 0, 2, fail_bin)
    ) %>%
    filter(
      futime > 0
    )
  
  rsample::initial_split(df, prop = 0.8, strata = fail_cr)
}

format_df_alldem <- function(df = df){
  df <- df %>%
    mutate(
      exit_date = exit_date_alldem,
      futime = futime_alldem,
      fail_bin = alldem_bin,
      fail_cr = ifelse(death == TRUE & fail_bin == 0, 2, fail_bin)
    ) %>%
    filter(
      futime > 0
    )
  
  rsample::initial_split(df, prop = 0.8, strata = fail_cr)
}