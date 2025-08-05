library(tidyverse)
library(lubridate)


# run helper script before loading data
source("~/RECOVER_SysBio/Utilities/helper_script.R")
env <- get_env_list("adult")
# all_forms <- env$formds_list() takes a long time, if not needed don't run

# Load the most recent data
files <- list.files("../project-files/Phenotype/")

pheno <- read.csv(paste0("../project-files/Phenotype/", files[length(files)]))

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

sd_cov <- pheno |> select(sd_income, sd_food, sd_homeless, sd_disability,
                          sd_docvisit, sd_unemploy, sd_medicaid, sd_uninsured,
                          sd_lostinsur, sd_moneyshort, sd_skipcare)


library(DescTools)
library(corrplot)

# Select only factor columns
n <- ncol(sd_cov)

# Initialize empty matrix
cramer_matrix <- matrix(NA, n, n, 
                        dimnames = list(names(sd_cov), names(sd_cov)))

# Compute pairwise Cramer’s V
for(i in 1:n) {
  for(j in 1:n) {
    cramer_matrix[i,j] <- CramerV(table(sd_cov[[i]], sd_cov[[j]]))
  }
}

head(cramer_matrix)

corrplot(cramer_matrix, is.corr = TRUE, method = "color",
         addCoef.col = "black", tl.col = "black", tl.srt = 45)

### Outcomes
# Symptom count
hist(pheno$symp_ct)
mod_symp_ct <- glm(symp_ct ~ age_enroll + as.factor(biosex_f) + bmi +
                     as.factor(alco_post) + as.factor(tobac_post),
                   data = pheno, family = poisson(link = "log"))

summary(mod_symp_ct)

# PASC score
hist(pheno$pasc_score_2024)
mod_pasc_score <- glm(pasc_score_2024 ~ age_enroll + as.factor(biosex_f) + bmi +
                        as.factor(alco_post) + as.factor(tobac_post),
                      family = poisson(link = "log"),
                      data = pheno)

summary(mod_pasc_score)

# Case (LC)
pheno$pasc_score_yn_2024 <- as.factor(pheno$pasc_score_yn_2024)
mod_lc_diag <- glm(pasc_score_yn_2024 ~ age_enroll + as.factor(biosex_f) + bmi +
                     as.factor(alco_post) + as.factor(tobac_post)
                   ,data = pheno, family = binomial(link = "logit"))

summary(mod_lc_diag)

# Socio economic status
# Hospitalized during acute phase

