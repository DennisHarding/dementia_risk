install.packages('readr')
sinstall.packages('tidyverse')
library(readr)
library(tidyverse)
library(data.table)



file_name1 <- file.choose()
file_name2 <- file.choose()

ukb_data <- fread(file_name1)
extended_data <- fread(file_name2)

extended_data_sel <- extended_data %>% select(ID, phylo_score)

censor_date_new <- max(ukb_data$visit_date_i3, na.rm = TRUE)

ukb_data <- tibble(ukb_data)

df <- ukb_data %>% 
  select(ID,
         entry_date,
         age_baseline,
         
         death_date,
         ad_hes_df,
         vd_hes_df,
         nonad_hes_df,
         
         gene_apoe,
         sex,
         British_white,
         
         alldm_hes_df,
         t2dm_hes_df,
         giga_ih_hes_df,
         is_hes_df,
         giga_stroke_hes_df,
         hf_hes_df,
         ihd_hes_df,
         ht_hes_df,
         wmh_total_i2,
         wmh_total_i3,
         chol,
         hdl,
         dbp,
         sbp,
         
         smok_ever,
         smok_stat,
         alc_freq,
         alc_gram,
         alcohol_unit_week,
         pa_IPAQ_group,
         edu_yrs
  ) %>%
  left_join(extended_data_sel %>% select(ID, phylo_score), by = "ID"
  ) 

write_tsv(df, "custom_data_X.tsv")
