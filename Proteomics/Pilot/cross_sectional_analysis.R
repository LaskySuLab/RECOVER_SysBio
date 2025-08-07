rm(list=ls())
library(dplyr)


source("/udd/spaap/Projects/RECOVER/code/RECOVER_SysBio/Utilities/aux_functions.R")

#define phenotype to import
phenos <- c("prote_topmed_id", "symp_ct", "pasc_score_2024", "pasc_score_yn_2024", 
            "age_enroll", "biosex_f", "bmi", "alco_pre", "alco_post", "tobac_pre", 
            "tobac_post", "timepoint", "participant_type", "part_id")

#read merged data
data <- merge_prote_pheno(phenos)

data <- data %>% 
  rename(age = age_enroll, 
  sex = biosex_f) %>%
  mutate(LC = ifelse(pasc_score_yn_2024 == "Yes", 1, 0),
         timepoint = factor_tp(timepoint)) %>%
  filter(!is.na(timepoint))

#initialize containers
LC_res <- list()
PS_res <- list()
SC_res <- list()

for (tp in tps) {
  message("Processing timepoint: ", tp)
  
  #subset data for the timepoint
  data_tp <- data %>%
    filter(timepoint == tp)
  
  n_samples <- nrow(data_tp)
  n_cases <- sum(data_tp$LC, na.rm = TRUE)
  
  #initialize containers
  resLC_t <- list()
  resPS_t <- list()
  resSC_t <- list()
  
  for (protein in prote_id) {
    
    #subset data for protein
    df <- data_tp %>%
      select(all_of(c(protein, "LC", "pasc_score_2024", "symp_ct", "age", "sex")))
    
    #Long Covid y/n (log reg)
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
    
    #PASC score (lin reg)
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
    
    #Symptom count (lin reg)
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
  
  #gather results
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

LC_res <- bind_rows(LC_res)
PS_res <- bind_rows(PS_res)
SC_res <- bind_rows(SC_res)


head(LC_res)
min(LC_res$fdr)
min(PS_res$fdr)
min(SC_res$fdr)

write.csv(LC_res,"/udd/spaap/Projects/RECOVER/data/proteomics/pilot/pilot_LC_logreg_age_sex_adj.csv", row.names = F)
write.csv(PS_res,"/udd/spaap/Projects/RECOVER/data/proteomics/pilot/pilot_PS_linreg_age_sex_adj.csv", row.names = F)
write.csv(SC_res,"/udd/spaap/Projects/RECOVER/data/proteomics/pilot/pilot_SC_linreg_age_sex_adj.csv", row.names = F)
