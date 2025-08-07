library(dplyr)

#--factors
tps <- c("pre_inf","0m","3m","6m","1y","2y","3y")
factor_tp <- function(data){
  return(factor(data, levels = tps))
}

LC_cat <- c("likely LC",    "possible LC", "no LC", "censored")
factor_LC <- function(data){
  return(factor(data , levels = names(LC_cat)))
}

#--colors
#significance colors
sig_colors <- c("ns" = "gray70", "pval<0.05" = "#E69F00", "fdr<0.05" = "#6A0DAD")
#long covid categories
colors_LC <- setNames(c("#D36C6C", "#4A90E2", "#009E73", "white"), LC_cat)


#--plots
theme <- theme(
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 12),
  strip.text = element_text(size = 12),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12)
)

#--read data
read_prote <- function(){
#read QCd proteomics data
  prote <- read.csv("/proj/rerecs/rerec00/data/protein/RECOVER_TOPMed_2025/results/proteomics_qc/20250715_Pilot_2000/20250805/protein_20250805.csv")
  #protein ID
  prote_id <- colnames(prote)[grepl("OID", colnames(prote))]
  #select only protein and tomped id columns
  prote_s <- prote %>%
    select(TOPMed.ID, all_of(prote_id))
  return(prote_s)
}

read_pheno <-  function(phenos){
  #read phenotype
    pheno_all <- read.csv("/proj/rerecs/rerec00/data/phenotype/data/protein_pilot/phenotype_adult_prote_pilot_masked20250805.csv") 
  #subset phenotype columns
  pheno_s <- pheno_all %>%
    select(any_of(unique(c("prote_topmed_id", phenos))))  %>%
    rename(TOPMed.ID = prote_topmed_id)
  return(pheno_s)
}

merge_prote_pheno <- function(phenos){
  pheno_s <- read_pheno(phenos)
  data <- read_prote() %>%
    left_join(pheno_s, by = "TOPMed.ID") %>%
    filter(TOPMed.ID %in% pheno_s$TOPMed.ID)
}


