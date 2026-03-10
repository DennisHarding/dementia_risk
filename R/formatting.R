#> ------------------------------------
#> FUNCTIONS FOR DATA FORMATTING
#> ------------------------------------

format_df <- function(data = data){
  # Update sensor date
  censor_date_new <- as.IDate("2026-01-01")

  df <- data %>%
    # Shorten variable names
    rename(
      ad_hes = ad_hes_df,
      vd_hes = vd_hes_df,
      nonad_hes = nonad_hes_df,   # vascular-related dementia
      prs = phylo_score
    ) %>%
    mutate(
      
      # Define exit dates (earliest of event, death, or censoring)
      exit_date_ad = pmin(death_date, ad_hes, censor_date_new, na.rm=TRUE),
      exit_date_vd = pmin(death_date, vd_hes, censor_date_new, na.rm=TRUE),
      exit_date_nonad = pmin(death_date, nonad_hes, censor_date_new, na.rm=TRUE),
      exit_date_alldem = pmin(death_date, ad_hes, vd_hes, nonad_hes, censor_date_new, na.rm=TRUE),
      
      # Follow-up time (years)
      futime_ad = as.numeric(difftime(exit_date_ad, entry_date, unit = "days")) / 365.25,
      futime_vd = as.numeric(difftime(exit_date_vd, entry_date, unit = "days")) / 365.25,
      futime_nonad = as.numeric(difftime(exit_date_nonad, entry_date, unit = "days")) / 365.25,
      futime_alldem = as.numeric(difftime(exit_date_alldem, entry_date, unit = "days")) / 365.25,
      
      # Age at event
      age_ad = age_baseline + as.numeric(difftime(ad_hes, entry_date, unit = "days")) / 365.25,
      age_vd = age_baseline + as.numeric(difftime(vd_hes, entry_date, unit = "days")) / 365.25, 
      age_nonad = age_baseline + as.numeric(difftime(nonad_hes, entry_date, unit = "days")) / 365.25, 
      
      # Event indicators
      ad_bin = as.numeric(!is.na(ad_hes)),
      vd_bin = as.numeric(!is.na(vd_hes)),
      nonad_bin = as.numeric(!is.na(nonad_hes)),
      alldem_bin = as.numeric(!is.na(coalesce(age_ad, age_vd, age_nonad))),
      
      # Education variables
      edu = case_when(
        edu_yrs == -3 ~ NA_real_,   # missing code
        edu_yrs >= 9 ~ 1,
        edu_yrs < 9 ~ 0,
        .default = NA_real_
      ),
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
      ),
      
      # Alcohol frequency grouped
      alc = case_when(
        alc_freq %in% c(1,2,3) ~ 0,
        alc_freq %in% c(4,5,6) ~ 1,
        .default = NA_real_
      ),
      
      # Earliest stroke record
      stroke_first = pmin(giga_stroke_hes_df, giga_ih_hes_df, is_hes_df, na.rm=TRUE),
      
      # Baseline disease indicators (present before entry)
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
      
      # Derived lipid measure
      non_hdl = chol - hdl,
      
      # Set APOE reference genotype
      gene_apoe = relevel(as.factor(gene_apoe), ref = "e33")
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
  df_ad
}

# Vascular Dementia
format_df_vd <- function(df = df){
  df %>%
    mutate(
      exit_date = exit_date_vd,
      futime = futime_vd,
      fail_bin = vd_bin,
      fail_age = age_vd
    ) %>%
    filter(
      futime > 0
    )
}

# Vascular-Related Dementia
format_df_nonad <- function(df = df){
  df %>%
    mutate(
      exit_date = exit_date_nonad,
      futime = futime_nonad,
      fail_bin = nonad_bin,
      fail_age = age_nonad
    ) %>%
    filter(
      futime > 0
    )
}

format_df_alldem <- function(df = df){
  df %>%
    mutate(
      exit_date = exit_date_alldem,
      futime = futime_alldem,
      fail_bin = alldem_bin,
    ) %>%
    filter(
      futime > 0
    )
}