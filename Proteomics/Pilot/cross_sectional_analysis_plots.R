library(ggridges)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggrepel)

#load protein info
prote_info <- read.csv("/proj/rerecs/rerec00/data/protein/RECOVER_TOPMed_2025/results/proteomics_qc/20250715_Pilot_2000/20250805/protein_info_20250805.csv")

#load and annotate significance
standardize_res <- function(res, outcome_label) {
  res %>%
    rename(OlinkID = protein) %>%
    left_join(prote_info %>% select(OlinkID, Assay), by = "OlinkID") %>%
    mutate(
      timepoint = factor_tp(timepoint),
      outcome = outcome_label,
      sig = case_when(
        fdr < 0.05 ~ "fdr<0.05",
        p_value < 0.05 ~ "pval<0.05",
        TRUE ~ "ns"
      ),
      log10p = -log10(p_value),
      label = ifelse(log10p > 3.5, Assay, NA)
    )
}

#load results
LC_res <- read.csv("/udd/spaap/Projects/RECOVER/data/proteomics/pilot/pilot_LC_logreg_age_sex_adj.csv") %>%
  standardize_res("LC")

PS_res <- read.csv("/udd/spaap/Projects/RECOVER/data/proteomics/pilot/pilot_PS_linreg_age_sex_adj.csv") %>%
  standardize_res("PS")

SC_res <- read.csv("/udd/spaap/Projects/RECOVER/data/proteomics/pilot/pilot_SC_linreg_age_sex_adj.csv") %>%
  standardize_res("SC")

all_res <- bind_rows(LC_res, PS_res, SC_res) %>%
  filter(timepoint != "pre_inf")

# Plot
png("/udd/spaap/Projects/RECOVER/figures/proteomics/pilot/pilot_LC_PS_SC_volcano_by_tp.png",
    type = "cairo", height = 900, width = 1400, res=100)

ggplot(all_res, aes(x = estimate, y = log10p, color = sig)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), color = "gray23", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray23", linetype = "dotted") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 100) +
  scale_color_manual(values = sig_colors, name = "Significance") +
  facet_grid(outcome ~ timepoint, scales = "free") +
  labs(y = expression(-log[10](p)), x = "Estimate") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "none"
  )

dev.off()


png("/udd/spaap/Projects/RECOVER/figures/proteomics/pilot/pilot_LC_logreg_age_sex_adj_pvalsdist.png", type="cairo", height=800, width = 1000)
ggplot(LC_res, aes(x = -log10(p_value), y = timepoint, fill = timepoint)) +
  geom_density_ridges() +
  geom_vline(xintercept = -log10(0.05))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off()
