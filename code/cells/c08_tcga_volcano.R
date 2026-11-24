# ==============================================================================
# Cell 08 - TCGA-CESC: differential expression, Tumour vs Normal, and volcano plot
# (manuscript Fig. 2; file Fig6_volcano_TCGA_CESC)
# Model: limma-trend on log2(TPM + 1) (Ritchie et al. 2015, Nucleic Acids Res; Law et al. 2014, Genome Biol),
# with robust empirical-Bayes moderation. Only 3 normal samples exist, so the variance moderation matters;
# results are descriptive for the normal arm. Unbalanced design: 304 tumours vs 3 normals.
# Needs in the session: expr_t, tcga_group (from cell 06) and the helpers of cell 01.
# ==============================================================================
keep <- rowSums(expr_t > 1) >= 3                      # TPM > 1 in at least 3 samples
log_num("tcga_deg_n_tested", sum(keep))
m8 <- expr_t[keep, ]
design8 <- model.matrix(~ 0 + tcga_group); colnames(design8) <- c("Normal", "Tumor")
fit8 <- lmFit(m8, design8)
fit8 <- eBayes(contrasts.fit(fit8, makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design8)),
               trend = TRUE, robust = TRUE)
deg_t <- topTable(fit8, coef = "Tumor_vs_Normal", number = Inf, adjust.method = "BH", sort.by = "none") %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::mutate(class = dplyr::case_when(adj.P.Val < FDR_CUT & logFC >  LFC_CUT ~ "Up in tumour",
                                         adj.P.Val < FDR_CUT & logFC < -LFC_CUT ~ "Down in tumour",
                                         TRUE ~ "Not significant"))
write.csv(deg_t, "results/tables/TableS6_TCGA_limma_Tumor_vs_Normal.csv", row.names = FALSE)
n_up_t <- sum(deg_t$class == "Up in tumour"); n_dn_t <- sum(deg_t$class == "Down in tumour")
log_num("tcga_deg_n_up", n_up_t); log_num("tcga_deg_n_down", n_dn_t)
log_num("tcga_deg_logFC_min", round(min(deg_t$logFC), 2)); log_num("tcga_deg_logFC_max", round(max(deg_t$logFC), 2))
log_num("tcga_deg_max_neglog10FDR", round(max(-log10(deg_t$adj.P.Val)), 2))
cat(sprintf("TCGA tumour vs normal: %d up, %d down (FDR<%.2f, |log2FC|>%d), %d genes tested\n",
            n_up_t, n_dn_t, FDR_CUT, LFC_CUT, nrow(deg_t)))
# Concordance with the cervical progression data (GSE63514, Cancer vs Normal), genes in common
if (exists("deg")) {
  cm <- merge(deg[, c("gene", "logFC", "class")], deg_t[, c("gene", "logFC", "class")], by = "gene", suffixes = c("_gse", "_tcga"))
  log_num("tcga_vs_gse_n_common_genes", nrow(cm))
  log_num("tcga_vs_gse_logFC_spearman", round(cor(cm$logFC_gse, cm$logFC_tcga, method = "spearman"), 3))
  both <- cm$class_gse != "Not significant" & cm$class_tcga != "Not significant"
  log_num("tcga_vs_gse_n_both_sig", sum(both))
  log_num("tcga_vs_gse_n_both_sig_same_direction", sum(both & sign(cm$logFC_gse) == sign(cm$logFC_tcga)))
}
lab_t <- bind_rows(deg_t %>% filter(class == "Up in tumour")   %>% slice_max(logFC, n = 10),
                   deg_t %>% filter(class == "Down in tumour") %>% slice_min(logFC, n = 10))
write.csv(lab_t, "results/tables/Fig6_labelled_genes.csv", row.names = FALSE)
cols_t <- c("Up in tumour" = "#D62728", "Down in tumour" = "#1F77B4", "Not significant" = "grey70")
p6 <- ggplot(deg_t, aes(logFC, -log10(adj.P.Val), colour = class)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = lab_t, aes(label = gene), size = 2.8, colour = "black", fontface = "italic",
                  max.overlaps = Inf, box.padding = 0.35, segment.size = 0.25, min.segment.length = 0) +
  scale_colour_manual(values = cols_t, name = NULL) +
  labs(x = expression(log[2]~"fold change (Tumour / Normal)"), y = expression(-log[10]~"FDR")) +
  theme_paper + theme(legend.position = "bottom")
save_fig(p6, "Fig6_volcano_TCGA_CESC", 7, 5.5)
print(p6)
