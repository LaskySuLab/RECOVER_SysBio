library(tidyverse)
library(lubridate)

file <- list.files("/sbgenomics/project-files/Phenotype/")

pheno <- read.csv("/sbgenomics/project-files/Phenotype/phenotype_adult_masked_20250728.csv")
