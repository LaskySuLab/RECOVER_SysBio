rm(list=ls())

source("/udd/spaap/Projects/RECOVER/code/RECOVER_SysBio/Utilities/aux_functions.R")

phenos <- c("prote_topmed_id", "symp_ct", "pasc_score_2024", "pasc_score_yn_2024", 
            "age_enroll", "biosex_f", "bmi", "alco_pre", "alco_post", "tobac_pre", 
            "tobac_post", "timepoint", "participant_type", "part_id")

data <- merge_prote_pheno(phenos)

sex_prote <- c("SHBG",	"LEP","ESR1","KDM5D","IGFBP1", "IGFBP2")
sex_prote_info <- read.csv("/proj/rerecs/rerec00/data/protein/RECOVER_TOPMed_2025/results/proteomics_qc/20250715_Pilot_2000/20250805/protein_info_20250805.csv") %>%
  filter(Assay %in% sex_prote)
data_s <- data %>%
  select(any_of(c("biosex_f",sex_prote_info$OlinkID)))
 
data_long <- data_s %>%
  pivot_longer(
    cols = -biosex_f,
    names_to = "OlinkID",
    values_to = "NPX"
  ) %>%
  left_join(sex_prote_info, by = "OlinkID") %>%
  mutate(Assay = fct_relevel(Assay, sex_prote))  %>%
  filter(!is.na(biosex_f))

png("/udd/spaap/Projects/RECOVER/figures/proteomics/pilot/pilot_sex_protes.png", type="cairo", height=800, width = 800)
ggplot(data_long, aes(x = biosex_f, y = NPX, fill = biosex_f)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~Assay, scales = "free_y", nrow = 2) +
  labs(x = "Sex", y = "NPX (Normalized Protein Expression)", title = "Protein Levels by Sex") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
dev.off()
