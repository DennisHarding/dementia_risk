plot_forest <- function(type, data){

    names <- c("alldm", 
               "smok_ever", 
               "edu", 
               "alc",
               "is_ih",
               "sbp10",
               "dbp10",
               "ht",
               "ihd"
    )
    
    legends <- list("alldm" = "Diabetes", 
                    "smok_ever" = "Smoking", 
                    "edu" = "Education", 
                    "alc" = "Alcohol intake",
                    "is_ih" = "is_ih stroke",
                    "sbp10" = "Systolic BP 10",
                    "dbp10" = "Diastolic BP 10",
                    "ht" = "Hypertension",
                    "ihd" = "Ischemic Heart Disease"
    )
    chain <- str_flatten(names, collapse = " + ")
    term <- c()
    labels_term <- c()
    conf_min_uni <- c()
    conf_max_uni <- c()
    estimate_uni <- c()
    conf_min_uni_gen <- c()
    conf_max_uni_gen <- c()
    estimate_uni_gen <- c()
    conf_min_all <- c()
    conf_max_all <- c()
    estimate_all <- c()
    conf_min_all_gen <- c()
    conf_max_all_gen <- c()
    estimate_all_gen <- c()
    
    formula_all <- paste0("Surv(futime, event = fail_bin) ~ sex + I(age_baseline^2) + ", chain)
    cox_all <- suppressWarnings(coxph(as.formula(formula_all), data = data))
    model_all <- tidy(cox_all, conf.int = TRUE, exponentiate = TRUE)
    formula_all_gen <- paste0("Surv(futime, event = fail_bin) ~ sex + I(age_baseline^2) + prs_fac + gene_apoe + ", chain)
    cox_all_gen <- suppressWarnings(coxph(as.formula(formula_all_gen), data = data))
    model_all_gen <- tidy(cox_all_gen, conf.int = TRUE, exponentiate = TRUE)
    
    
    for (i in names){
      term <- c(term, i)
      labels_term <- c(labels_term, legends[[i]])
      conf_min_all <- c(conf_min_all,as.numeric(model_all[str_detect(model_all$term, i),]$conf.low))
      conf_max_all <- c(conf_max_all,as.numeric(model_all[str_detect(model_all$term, i),]$conf.high))
      estimate_all <- c(estimate_all,as.numeric(model_all[str_detect(model_all$term, i),]$estimate))
      conf_min_all_gen <- c(conf_min_all_gen,as.numeric(model_all_gen[str_detect(model_all_gen$term, i),]$conf.low))
      conf_max_all_gen <- c(conf_max_all_gen,as.numeric(model_all_gen[str_detect(model_all_gen$term, i),]$conf.high))
      estimate_all_gen <- c(estimate_all_gen,as.numeric(model_all_gen[str_detect(model_all_gen$term, i),]$estimate))
      formula_uni <- paste0("Surv(futime, event = fail_bin) ~ sex + I(age_baseline^2) + ", i)
      cox_uni <- suppressWarnings(coxph(as.formula(formula_uni), data = data))
      model_uni <- tidy(cox_uni, conf.int = TRUE, exponentiate = TRUE)
      conf_min_uni <- c(conf_min_uni,as.numeric(model_uni[str_detect(model_uni$term, i),]$conf.low))
      conf_max_uni <- c(conf_max_uni,as.numeric(model_uni[str_detect(model_uni$term, i),]$conf.high))
      estimate_uni <- c(estimate_uni,as.numeric(model_uni[str_detect(model_uni$term, i),]$estimate))
      formula_uni_gen <- paste0("Surv(futime, event = fail_bin) ~ sex + I(age_baseline^2) + prs_fac + gene_apoe + ", i)
      cox_uni_gen <- suppressWarnings(coxph(as.formula(formula_uni_gen), data = data))
      model_uni_gen <- tidy(cox_uni_gen, conf.int = TRUE, exponentiate = TRUE)
      conf_min_uni_gen <- c(conf_min_uni_gen,as.numeric(model_uni_gen[str_detect(model_uni_gen$term, i),]$conf.low))
      conf_max_uni_gen <- c(conf_max_uni_gen,as.numeric(model_uni_gen[str_detect(model_uni_gen$term, i),]$conf.high))
      estimate_uni_gen <- c(estimate_uni_gen,as.numeric(model_uni_gen[str_detect(model_uni_gen$term, i),]$estimate))
    }
    
    prs_rows <- model_all_gen %>%
      filter(str_detect(term, "prs_fac")) %>%
      mutate(term = paste0("prs", substr(term, nchar(term), nchar(term)), " vs prs3"))
    
    apoe_rows <- model_all_gen %>%
      filter(str_detect(term, "gene_apoe")) %>%
      mutate(term = paste0(substr(term, nchar(term) - 2, nchar(term)), " vs e33"))
    
    prs_df <- prs_rows %>%
      transmute(
        term = term,
        conf_min = conf.low,
        conf_max = conf.high,
        estimate = estimate,
        model = "all_gen"
      )
    
    apoe_df <- apoe_rows %>%
      transmute(
        term = term,
        conf_min = conf.low,
        conf_max = conf.high,
        estimate = estimate,
        model = "all_gen"
      )
    
    prs_header <- tibble(
      term = "Polygenetic Risk Score",
      conf = "",
      model = "",
      estimate = NA,
      conf_min = NA,
      conf_max = NA,
      is_summary = TRUE
    )
    
    apoe_header <- tibble(
      term = "ApoE variants",
      conf = "",
      model = "",
      estimate = NA,
      conf_min = NA,
      conf_max = NA,
      is_summary = TRUE
    )
    
    data_uni <- data.frame(
      term = labels_term, 
      conf_min = conf_min_uni, 
      conf_max = conf_max_uni, 
      estimate = estimate_uni) %>% 
      mutate(model = "uni")
    data_uni_gen <- data.frame(
      term = labels_term, 
      conf_min = conf_min_uni_gen, 
      conf_max = conf_max_uni_gen, 
      estimate = estimate_uni_gen) %>% 
      mutate(model = "uni_gen")
    data_all <- data.frame(
      term = labels_term, 
      conf_min = conf_min_all, 
      conf_max = conf_max_all, 
      estimate = estimate_all) %>% 
      mutate(model = "all")
    data_all_gen <- data.frame(
      term = labels_term, 
      conf_min = conf_min_all_gen, 
      conf_max = conf_max_all_gen, 
      estimate = estimate_all_gen) %>% 
      mutate(model = "all_gen")
    df_forest <- bind_rows(data_uni, data_uni_gen, data_all, data_all_gen) 
    
    df_plot <- df_forest %>%
      arrange(match(term, unlist(legends[names])), match(model, c("uni", "uni_gen", "all", "all_gen"))) %>%
      mutate(
        conf = paste0(sprintf("%.2f", estimate), " (", sprintf("%.2f", conf_min), ", ", sprintf("%.2f", conf_max), ")"),
        term = fct_inorder(term)
      ) %>%
      group_split(term) %>%
      lapply(function(df_term) {
        term_name <- unique(df_term$term)
        header <- tibble(term = term_name, 
                         conf = "", 
                         model = "", 
                         estimate = NA, 
                         conf_min = NA, 
                         conf_max = NA, 
                         is_summary = TRUE)
        rows   <- mutate(df_term, is_summary = FALSE)
        spacer <- tibble(term = "", 
                         conf = "", 
                         model = "", 
                         estimate = NA, 
                         conf_min = NA, 
                         conf_max = NA, 
                         is_summary = FALSE)
        bind_rows(header, rows, spacer)
      }) %>%
      bind_rows() %>%
      mutate(
        term = case_when(
          (!is_summary) & (model == "uni") ~ "Sex + Age",
          (!is_summary) & (model == "uni_gen") ~ "Sex + Age w/ genetics",
          (!is_summary) & (model == "all") ~ "All covariates",
          (!is_summary) & (model == "all_gen") ~ "All covariates w/ genetics",
          is_summary ~ term,
          .default = NA
        ),
        term = fct_inorder(term)
      )
    
    prs_rows <- bind_rows(prs_df) %>%
      mutate(
        conf = paste0(sprintf("%.2f", estimate), 
                      " (", sprintf("%.2f", conf_min), 
                      ", ", sprintf("%.2f", conf_max), ")"),
        is_summary = FALSE
      )
    apoe_rows <- bind_rows(apoe_df) %>%
      mutate(
        conf = paste0(sprintf("%.2f", estimate), 
                      " (", sprintf("%.2f", conf_min), 
                      ", ", sprintf("%.2f", conf_max), ")"),
        is_summary = FALSE
      )
    spacer <- tibble(term = "", 
                     conf = "", 
                     model = "", 
                     estimate = NA, 
                     conf_min = NA, 
                     conf_max = NA, 
                     is_summary = FALSE)
    
    df_plot <- bind_rows(df_plot, apoe_header, apoe_rows, spacer, prs_header, prs_rows)
    
    risk_factors_header <- tibble(
      term = paste("Risk factor -", type), 
      conf = "HR (95% CI)", 
      model = "", 
      estimate = NA, 
      conf_min = NA, 
      conf_max = NA, 
      is_summary = TRUE
    )
    df_plot <- bind_rows(risk_factors_header, df_plot)
    
    forest <- forestplot(
      labeltext = cbind(df_plot$term, df_plot$conf),
      mean      = df_plot$estimate,
      lower     = df_plot$conf_min,
      upper     = df_plot$conf_max,
      is.summary = df_plot$is_summary,
      zero = TRUE,
      lwd.zero = 1,
      boxsize = 0.3,
      xlab = "HR (95% CI)",
      xlog = TRUE,
      xticks = c(0.35, 0.5, 1, 2, 4, 6, 8),
      ci.vertices = TRUE,
      fn.ci_norm = fpDrawCircleCI,
      lwd.ci = 2,
      line.margin = 1,
      col = fpColors(box = "black", line = "black", zero = "black"),
      txt_gp = fpTxtGp(cex = 2, xlab = gpar(cex = 2, fontface = "bold"), ticks = gpar(cex = 2)),
      mar = unit(c(0.5, 1, 0.5, 1), "cm")
    )
    png(paste0("figures/forestplot/forest_", type, ".png"), width=18, height=22.25, units="in", res=300)
    print(forest)
    forest
    dev.off()
  }
