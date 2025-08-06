library(ggridges)
library(ggplot2)
library(ggrepel)

prote_info <- read.csv("/proj/rerecs/rerec00/data/protein/RECOVER_TOPMed_2025/results/proteomics_qc/20250715_Pilot_2000/20250805/protein_info_20250805.csv")

LC_res <- read.csv("/udd/spaap/Projects/RECOVER/data/proteomics/pilot/pilot_LC_logreg_age_sex_adj.csv")

LC_res <- LC_res %>%
  rename(OlinkID = protein)%>%
  left_join(prote_info %>% select(OlinkID ,Assay), by = "OlinkID")%>%
  mutate(timepoint = factor_tp(timepoint),
         sig = ifelse(p_value < 0.05, "pval<0.05", "ns"),
         sig = ifelse(fdr < 0.1, "fdr<0.1", sig),
         label = ifelse(-log10(p_value) > 3.5, Assay, NA))
head(LC_res)

png("/udd/spaap/Projects/RECOVER/figures/proteomics/pilot/pilot_LC_logreg_age_sex_adj_pvalsdist.png", type="cairo", height=800, width = 800)
ggplot(LC_res, aes(x = -log10(p_value), y = timepoint, fill = timepoint)) +
  geom_density_ridges() +
  geom_vline(xintercept = -log10(0.05))+
  theme_ridges() + 
  theme(legend.position = "none")
dev.off()


png("/udd/spaap/Projects/RECOVER/figures/proteomics/pilot/pilot_LC_logreg_age_sex_adj_volcano_tp.png", type="cairo", height=300, width = 1000)
ggplot(LC_res %>% filter(timepoint != "pre_inf"), 
       aes(y = -log10(p_value), x = estimate, color = sig)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), color = "gray23", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray23", linetype = "dotted") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 100) +
  facet_grid(~timepoint, scales = "free") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )
dev.off()

