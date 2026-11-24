# ==============================================================================
# Cell 03 - Differential expression (limma) and volcano plot (Fig. 2)
# ==============================================================================
run_limma <- function(mat, group, covariate = NULL) {
  design <- if (is.null(covariate)) model.matrix(~ 0 + group) else model.matrix(~ 0 + group + covariate)
  colnames(design) <- make.names(sub("^group", "", colnames(design)))
  fit <- lmFit(mat, design)
  cm  <- makeContrasts(Cancer_vs_Normal = Cancer - Normal, CIN1_vs_Normal = CIN1 - Normal,
                       CIN2_vs_Normal = CIN2 - Normal, CIN3_vs_Normal = CIN3 - Normal,
                       Cancer_vs_CIN3 = Cancer - CIN3, levels = design)
  eBayes(contrasts.fit(fit, cm))               # empirical-Bayes moderation
}
classify <- function(tab) dplyr::mutate(tab, class = dplyr::case_when(
  adj.P.Val < FDR_CUT & logFC >  LFC_CUT ~ "Up in cancer",
  adj.P.Val < FDR_CUT & logFC < -LFC_CUT ~ "Down in cancer", TRUE ~ "Not significant"))
get_tab <- function(fit, coef) topTable(fit, coef = coef, number = Inf, adjust.method = "BH", sort.by = "none") %>%
  tibble::rownames_to_column("gene")

fit  <- run_limma(expr, group)
deg  <- classify(get_tab(fit, "Cancer_vs_Normal"))
write.csv(deg, "results/tables/TableS2_GSE63514_limma_Cancer_vs_Normal.csv", row.names = FALSE)
for (cf in c("CIN1_vs_Normal", "CIN2_vs_Normal", "CIN3_vs_Normal", "Cancer_vs_CIN3"))
  write.csv(classify(get_tab(fit, cf)), sprintf("results/tables/TableS2_GSE63514_limma_%s.csv", cf), row.names = FALSE)

n_up <- sum(deg$class == "Up in cancer"); n_dn <- sum(deg$class == "Down in cancer")
log_num("deg_n_tested", nrow(deg)); log_num("deg_n_up", n_up); log_num("deg_n_down", n_dn)
log_num("deg_logFC_min", round(min(deg$logFC), 2)); log_num("deg_logFC_max", round(max(deg$logFC), 2))
log_num("deg_max_neglog10FDR", round(max(-log10(deg$adj.P.Val)), 2))
cat(sprintf("Cancer vs Normal: %d up, %d down (FDR<%.2f, |log2FC|>%d)\n", n_up, n_dn, FDR_CUT, LFC_CUT))

# ---- Sensitivity analyses (reported in the Supplement) ------------------------
# (a) with MT-/ribosomal transcripts retained
deg_with <- classify(get_tab(run_limma(expr_all, group), "Cancer_vs_Normal"))
sig_main <- deg$gene[deg$class != "Not significant"]; sig_with <- deg_with$gene[deg_with$class != "Not significant"]
log_num("sens_with_MTribo_n_deg", length(sig_with)); log_num("sens_with_MTribo_overlap", length(intersect(sig_main, sig_with)))
# (b) adjusting for dissection method (laser capture vs whole section)
fit_d <- run_limma(expr, group, covariate = dissection)
deg_d <- classify(get_tab(fit_d, "Cancer_vs_Normal"))
sig_d <- deg_d$gene[deg_d$class != "Not significant"]
log_num("sens_dissection_n_deg", length(sig_d)); log_num("sens_dissection_overlap", length(intersect(sig_main, sig_d)))
log_num("sens_dissection_logFC_cor", round(cor(deg$logFC, deg_d$logFC[match(deg$gene, deg_d$gene)]), 3))

# ---- Fig. 2: volcano with the 10 strongest up- and down-regulated genes labelled
lab <- bind_rows(deg %>% filter(class == "Up in cancer")   %>% slice_max(logFC, n = 10),
                 deg %>% filter(class == "Down in cancer") %>% slice_min(logFC, n = 10))
write.csv(lab, "results/tables/Fig2_labelled_genes.csv", row.names = FALSE)
cols <- c("Up in cancer" = "#D62728", "Down in cancer" = "#1F77B4", "Not significant" = "grey70")
p2 <- ggplot(deg, aes(logFC, -log10(adj.P.Val), colour = class)) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = lab, aes(label = gene), size = 2.8, colour = "black", fontface = "italic",
                  max.overlaps = Inf, box.padding = 0.35, segment.size = 0.25, min.segment.length = 0) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(x = expression(log[2]~"fold change (Cancer / Normal)"), y = expression(-log[10]~"FDR")) +
  theme_paper + theme(legend.position = "bottom")
save_fig(p2, "Fig2_volcano_GSE63514", 7, 5.5)
print(p2)
