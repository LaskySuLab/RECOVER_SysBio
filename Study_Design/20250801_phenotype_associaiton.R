library(tidyverse)
library(DescTools)
library(corrplot)
library(RColorBrewer)
library(data.table)

setwd("/proj/rerecs/rerec00/data/phenotype/data/")
options(bitmapType='cairo')
# run helper script before loading data on SEVEN BRIDGES
# source("~/RECOVER_SysBio/Utilities/helper_script.R")
# env <- get_env_list("adult")
# all_forms <- env$formds_list() takes a long time, if not needed don't run

# Load the most recent data
pheno <- read.csv("./phenotype_adult_masked_20250805.csv")

# Basic associations of outcomes and covariates

# social determinants

pheno <- pheno |>
  mutate(sd_income = as.factor(sd_income),
         sd_food = as.factor(sd_food),
         sd_homeless = as.factor(sd_homeless),
         sd_disability = as.factor(sd_disability),
         sd_docvisit = as.factor(sd_docvisit),
         sd_unemploy = as.factor(sd_unemploy),
         sd_medicaid = as.factor(sd_medicaid),
         sd_uninsured = as.factor(sd_uninsured),
         sd_lostinsur = as.factor(sd_lostinsur),
         sd_moneyshort = as.factor(sd_moneyshort),
         sd_skipcare = as.factor(sd_skipcare)) 

# select only one entry at one timepoint for one individual 
# these varialbe should all be the same withinin individual though
sd_cov <- pheno |> group_by(part_id) |>
  summarise(across(c(sd_income, sd_food, sd_homeless, sd_disability,
                     sd_docvisit, sd_unemploy, sd_medicaid, sd_uninsured,
                     sd_lostinsur, sd_moneyshort, sd_skipcare),
                   first), .groups = "drop")

# double check we have all the unique individual
length(unique(pheno$part_id)) == nrow(sd_cov) # TRUE

# Select only factor columns
n <- ncol(sd_cov[, -1])

# Initialize empty matrix
cramer_matrix <- matrix(NA, n, n, 
                        dimnames = list(names(sd_cov[, -1]), names(sd_cov[, -1])))

# Compute pairwise Cramer’s V
for(i in 1:n) {
  for(j in 1:n) {
    cramer_matrix[i,j] <- CramerV(table(sd_cov[, -1][[i]], sd_cov[, -1][[j]]))
  }
}

head(cramer_matrix)
col=rev(brewer.pal(n = 8, name = "RdBu"))
corrplot(cramer_matrix, is.corr = TRUE, method = "color",
         col = col,
         addCoef.col = "black", tl.col = "black", tl.srt = 45)

### Outcomes
# Symptom count
hist(pheno$symp_ct)
# mod_symp_ct <- glm(symp_ct ~ age_enroll + as.factor(biosex_f) + bmi +
#                      as.factor(alco_pre)+  as.factor(tobac_pre) +
#                      as.factor(alco_post) + as.factor(tobac_post),
#                    data = pheno, family = poisson(link = "log"))
# 
# summary(mod_symp_ct)

### single variate analysis 
cov_continuous <- c("age_enroll", "bmi")
cov_categorical <- c("alco_pre", "alco_post", "tobac_pre", "tobac_post",
                     "sd_income")

# covariate type either categorical or continuous
# return summary statistics for association with symptom count,
# PASC score, and PASC diagonsis (binary)
sing_var_association <- function(outcome_name, dat = pheno){ 
  results <- lapply(cov_continuous, function(var){
    formula <- as.formula(paste0(outcome_name," ~ ", var))
    if(outcome_name == "pasc_score_yn_2024") {
      res <- glm(formula, data = dat, family = binomial(link = "logit"))
    }else{
      res <- glm(formula, data = dat, family = poisson(link = "log"))
    }
    summary(res)
  })
  
  outcome <- lapply(results, function(z) c(z$coeff[2,1], z$coeff[2,2], z$coeff[2,3], z$coeff[2,4]))
  outcome1 <- as.data.table(do.call(rbind, outcome))
  colnames(outcome1) <- c("Beta","SE", "Z", "P")
  outcome1[, Outcome := outcome_name]
  outcome1[, Covariate := cov_continuous]
  return(outcome1)
}

sing_var_association_cat <- function(outcome_name, dat = pheno){ 
  results <- lapply(cov_categorical, function(var){
    # formula for ANOVA
    formula <- as.formula(paste0(outcome_name," ~ ", var))
    
    if(outcome_name == "pasc_score_yn_2024") {
      # Chi-sq test for binary outcome
      covariate <- as.factor(dat[[var]])
      pval <- chisq.test(dat[[outcome_name]], covariate)$p.value
    } else {
      # ANOVA for continuous outcome
      res <- aov(formula, data = dat)
      print(res)
      pval <- summary(res)[[1]][1, "Pr(>F)"]
    }
    
    # Return a row as list
    data.frame(
      Outcome = outcome_name,
      Covariate = var,
      P = pval
    )
  })
  
  # Combine list of data.frames into one data.frame
  results_df <- do.call(rbind, results)
  rownames(results_df) <- NULL
  
  return(results_df)
}


create_asso_pval_table <- function(dat){
  pheno$pasc_score_yn_2024 <- as.factor(pheno$pasc_score_yn_2024)
  asso_res_total <- c()
  for(outcome in c("symp_ct", "pasc_score_2024", "pasc_score_yn_2024")){
    asso_res_total <- rbind(asso_res_total, sing_var_association(outcome, dat))
  }
  
  asso_res_total_cat <- c()
  # pheno$sd_income <- as.factor(pheno$sd_income)
  for(outcome in c("symp_ct", "pasc_score_2024", "pasc_score_yn_2024")){
    asso_res_total_cat <- rbind(asso_res_total_cat, sing_var_association_cat(outcome, dat))
  }
  
  asso_res_total <- asso_res_total |> select(Outcome, Covariate, P)
  asso_res_total <- rbind(asso_res_total, asso_res_total_cat)
  asso_res_total_wide <- asso_res_total |> pivot_wider(names_from = Covariate,
                                                       values_from = P)
  return(as.data.frame(asso_res_total_wide))
}

time0 <- pheno |> filter(timepoint == "0m")
table_time0 <- create_asso_pval_table(time0)

time3m <-  pheno |> filter(timepoint == "3m")
table_time3m <- create_asso_pval_table(time3m)

time6m <-  pheno |> filter(timepoint == "6m")
table_time6m <- create_asso_pval_table(time6m)

time9m <-  pheno |> filter(timepoint == "9m")
table_time9m <- create_asso_pval_table(time9m)

time1y <-  pheno |> filter(timepoint == "1y")
table_time1y <- create_asso_pval_table(time1y)

time2y <-  pheno |> filter(timepoint == "2y")
table_time2y <- create_asso_pval_table(time2y)

time3y <-  pheno |> filter(timepoint == "3y")
table_time3y <- create_asso_pval_table(time3y)


