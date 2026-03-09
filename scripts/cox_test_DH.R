install.packages('readr', 'tidyverse')

library(survival)
library(readr)
library(dplyr)
library(ggplot2)
library(data.table)

system("dx download ./extended_data_v2.tsv")

odf <- fread("extended_data_v2.tsv")
odf <- as_tibble(odf)

summary(df)

df <- odf %>%
  select(gene_apoe, ad_bin, phylo_score, age_baseline, 
         smok_stat, ad_date_first, entry_date, ) %>% 
  mutate(
    apoe = gene_apoe,
    Cov = factor(smok_stat, levels = c("Never", "Previous", "Current")),
    age = age_baseline + 
      floor(
        as.numeric(difftime(
          ad_date_first, entry_date, unit = "days") / 365.25))) %>% 
  select(-gene_apoe)

model_cov <- coxph(Surv(time = age, event = ad_bin) ~ phylo_score + apoe
                   + Cov, data = df)


# Bringing in new data
start_age = c(50,60,70)
phylo_score = c(-1,-0.5,0,0.5,1)
apoe = sort(unique(df$apoe))
Cov = levels(df$Cov)
ml <- expand.grid(
  phylo_score = phylo_score,
  apoe = apoe,
  Cov = Cov)

# Cox Linear Risk Predictions
ml$pred <- predict(model_cov, newdata = ml, type="lp")
# Baseline approximation
base_surv <- survfit(model_cov)
base_df <- data.frame(time = base_surv$time, surv = base_surv$surv)
ml$base50 <- approx(base_df$time, base_df$surv, xout = 50)$y
ml$base60 <- approx(base_df$time, base_df$surv, xout = 60)$y
ml$base70 <- approx(base_df$time, base_df$surv, xout = 70)$y
ml$base80 <- approx(base_df$time, base_df$surv, xout = 80)$y


pred_cox <- data.frame(dis_pred_50 = 1 - (ml$base60/ml$base50)^(exp(ml$pred)),
                       dis_pred_60 = 1 - (ml$base70/ml$base60)^(exp(ml$pred)),
                       dis_pred_70 = 1 - (ml$base80/ml$base70)^(exp(ml$pred))) 
pred_cox <- pred_cox %>%
  mutate(col = case_when(
    dis_pred_70 < 0.01 ~ "<1",
    dis_pred_70 >= 0.01 & dis_pred_70 < 0.05 ~ "1",
    dis_pred_70 >= 0.05 & dis_pred_70 < 0.10 ~ "2",
    dis_pred_70 >= 0.10 & dis_pred_70 < 0.15 ~ "3",
    dis_pred_70 >= 0.15 & dis_pred_70 < 0.20 ~ "4",
    dis_pred_70 >= 0.20 & dis_pred_70 < 0.30 ~ "5",
    dis_pred_70 >= 0.30 ~ "6",
  ),
  text = ifelse(dis_pred_70 >= 0.20, "white", "black")
  )
cols <- c("<1"="#00994C", "1"="#66CC00","2"="#FFFF33","3"="#FF9933",
          "4"="#FF0000","5"="#CC0000", "6"="#994C00")
pred_cox <- cbind(ml, pred_cox)

p4 <- ggplot(pred_cox, aes(x = as.factor(apoe), y = as.factor(phylo_score),
                           fill = col)) +
  geom_tile() +
  geom_text(aes(label = lapply(dis_pred_70 * 100.0, as.integer), color = text),
            fontface = "bold") +
  scale_fill_manual(values = cols) +
  scale_color_identity() + 
  facet_grid(. ~ Cov) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed()

p4
ggsave("cox_model.png", plot = p4)


system("dx upload cox_test_DH.png --path dharding/")
