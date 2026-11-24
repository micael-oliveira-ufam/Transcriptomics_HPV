# ==============================================================================
# Cell 07 - TCGA-CESC: 18-gene panel, tumour vs normal (Wilcoxon, BH + Bonferroni)
# and hierarchical-clustering heatmap (Fig. 5)
# ==============================================================================
stopifnot(all(PANEL$gene %in% rownames(expr_t)))
tw <- lapply(PANEL$gene, function(g) {
  y <- expr_t[g, ]
  tibble(gene = g, median_normal = median(y[tcga_group == "Normal"]), median_tumor = median(y[tcga_group == "Tumor"]),
         delta_log2TPM = median(y[tcga_group == "Tumor"]) - median(y[tcga_group == "Normal"]),
         wilcox_p = wilcox.test(y ~ tcga_group, exact = TRUE)$p.value,          # exact permutation p (n = 3 vs 304)
         wilcox_p_asymptotic = wilcox.test(y ~ tcga_group, exact = FALSE)$p.value)  # normal approximation, for comparison
}) %>% bind_rows() %>% mutate(fdr = p.adjust(wilcox_p, "BH"), bonferroni = p.adjust(wilcox_p, "bonferroni")) %>%
  left_join(PANEL, by = "gene") %>% relocate(axis, .after = gene)
write.csv(tw, "results/tables/TableS5_TCGA_panel_tumor_vs_normal.csv", row.names = FALSE)
log_num("tcga_panel_n_bonferroni_sig", sum(tw$bonferroni < 0.05)); log_num("tcga_panel_n_bonferroni_sig_asymptotic", sum(p.adjust(tw$wilcox_p_asymptotic, "bonferroni") < 0.05)); log_num("tcga_panel_n_fdr_sig", sum(tw$fdr < 0.05))
n_t <- sum(tcga_group == "Tumor"); n_n <- sum(tcga_group == "Normal")
log_num("tcga_min_possible_exact_wilcoxon_p", signif(2 / choose(n_t + n_n, n_n), 3))
log_num("tcga_min_possible_bonferroni_p", signif(min(1, 18 * 2 / choose(n_t + n_n, n_n)), 3))
print(as.data.frame(tw %>% mutate(across(where(is.numeric), ~ signif(., 3))) %>% dplyr::select(gene, median_normal, median_tumor, wilcox_p, fdr, bonferroni)))

# Heatmap: row-wise Z-score truncated at +/-2, Euclidean distance, complete linkage
m <- expr_t[PANEL$gene, ]
z <- t(scale(t(m))); z[z > 2] <- 2; z[z < -2] <- -2
hc_cols <- hclust(dist(t(z), method = "euclidean"), method = "complete")   # same distance/linkage as pheatmap
cl <- cutree(hc_cols, k = 2)                                               # two sample clusters, named A and B
cl_lab <- factor(paste0("Cluster ", LETTERS[cl]))
log_num("tcga_heatmap_cluster_sizes", paste(table(cl_lab), collapse = "/"))
log_num("tcga_heatmap_normal_per_cluster", paste(table(cl_lab, tcga_group)[, "Normal"], collapse = "/"))
ann_col <- data.frame(Sample = tcga_group, Cluster = cl_lab, row.names = colnames(m))
ann_row <- data.frame(Axis = PANEL$axis, row.names = PANEL$gene)
ann_cols <- list(Sample = c(Normal = "#E69F00", Tumor = "#009E73"), Cluster = c("Cluster A" = "#CC79A7", "Cluster B" = "#56B4E9"))
write.csv(data.frame(sample = names(cl), cluster = cl, group = as.character(tcga_group)), "results/tables/Fig5_cluster_assignment.csv", row.names = FALSE)
for (ext in c("png", "tiff")) {
  f <- file.path("results/figures", paste0("Fig5_heatmap_TCGA_CESC.", ext))
  if (ext == "png") png(f, width = 3000, height = 1800, res = 300) else tiff(f, width = 3000, height = 1800, res = 300, compression = "lzw", type = "cairo")
  pheatmap(z, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100), clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean", clustering_method = "complete", annotation_col = ann_col, annotation_row = ann_row, annotation_colors = ann_cols,
           show_colnames = FALSE, fontsize_row = 9, treeheight_col = 30, treeheight_row = 30)
  dev.off()
}
