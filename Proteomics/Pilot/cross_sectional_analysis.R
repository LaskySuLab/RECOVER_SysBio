rm(list=ls())
library(dplyr)


source("/udd/spaap/Projects/RECOVER/code/RECOVER_SysBio/Utilities/aux_functions.R")

phenos <- c("prote_topmed_id", "symp_ct", "pasc_score_2024", "pasc_score_yn_2024", 
            "age_enroll", "biosex_f", "bmi", "alco_pre", "alco_post", "tobac_pre", 
            "tobac_post", "timepoint", "participant_type", "part_id")


data <- merge_prote_pheno(phenos)

data <- data %>% 
  rename(age = age_enroll, 
  sex = biosex_f) %>%
  mutate(LC = ifelse(pasc_score_yn_2024 == "Yes", 1, 0),
         timepoint = factor_tp(timepoint)) %>%
  filter(!is.na(timepoint))

LC_res <- list()
PS_res <- list()
SC_res <- list()

for (tp in tps) {
  message("Processing timepoint: ", tp)
  
  # Subset data for the timepoint
  data_tp <- data %>%
    filter(timepoint == tp)
  
  n_samples <- nrow(data_tp)
  n_cases <- sum(data_tp$LC, na.rm = TRUE)
  
  # Initialize containers
  resLC_t <- list()
  resPS_t <- list()
  resSC_t <- list()
  
  for (protein in prote_id) {
    
    df <- data_tp %>%
      select(all_of(c(protein, "LC", "pasc_score_2024", "symp_ct", "age", "sex")))
    
    #--- Logistic regression for LC
    form_LC <- as.formula(paste0("LC ~ ", protein, " + age + sex"))
    model_LC <- glm(form_LC, data = df, family = binomial)
    coef_LC <- summary(model_LC)$coefficients
    est_LC <- coef_LC[protein, "Estimate"]
    pval_LC <- coef_LC[protein, "Pr(>|z|)"]
    resLC_t[[protein]] <- data.frame(
      protein = protein,
      estimate = est_LC,
      p_value = pval_LC
    )
    
    #--- Linear regression for PS
    form_PS <- as.formula(paste0("pasc_score_2024 ~ ", protein, " + age + sex"))
    model_PS <- lm(form_PS, data = df)
    coef_PS <- summary(model_PS)$coefficients
    est_PS <- coef_PS[protein, "Estimate"]
    pval_PS <- coef_PS[protein, "Pr(>|t|)"]
    resPS_t[[protein]] <- data.frame(
      protein = protein,
      estimate = est_PS,
      p_value = pval_PS
    )
    
    #--- Linear regression for SC
    form_SC <- as.formula(paste0("symp_ct ~ ", protein, " + age + sex"))
    model_SC <- lm(form_SC, data = df)
    coef_SC <- summary(model_SC)$coefficients
    est_SC <- coef_SC[protein, "Estimate"]
    pval_SC <- coef_SC[protein, "Pr(>|t|)"]
    resSC_t[[protein]] <- data.frame(
      protein = protein,
      estimate = est_SC,
      p_value = pval_SC
    )
  }
  
  # Combine, add metadata, adjust p-values
  resLC_df <- bind_rows(resLC_t) %>%
    mutate(timepoint = tp, samples = n_samples, cases = n_cases,
           fdr = p.adjust(p_value, method = "fdr"))
  LC_res[[tp]] <- resLC_df
  
  resPS_df <- bind_rows(resPS_t) %>%
    mutate(timepoint = tp, samples = n_samples,
           fdr = p.adjust(p_value, method = "fdr"))
  PS_res[[tp]] <- resPS_df
  
  resSC_df <- bind_rows(resSC_t) %>%
    mutate(timepoint = tp, samples = n_samples,
           fdr = p.adjust(p_value, method = "fdr"))
  SC_res[[tp]] <- resSC_df
}

# Bind all timepoints
LC_res <- bind_rows(LC_res)
PS_res <- bind_rows(PS_res)
SC_res <- bind_rows(SC_res)





head(LC_res     )

LC_res_p <- LC_res %>%
  rename(OlinkID = protein)%>%
  left_join(prote_info %>% select(OlinkID ,Assay), by = "OlinkID")%>%
  mutate(timepoint = factor_tp(timepoint),
         sig = ifelse(p_value < 0.05, "pval<0.05", "ns"),
         sig = ifelse(fdr < 0.1, "fdr<0.1", sig),
         label = ifelse(-log10(p_value) > 3.5, Assay, NA))

