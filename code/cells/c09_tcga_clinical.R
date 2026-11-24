# ==============================================================================
# Cell 09 - TCGA-CESC clinical descriptors, annotated panel heatmap (Fig. 5) and
# panel expression by histology / FIGO stage / HPV genotype.
# Needs from the session: expr_t (log2 TPM, MT-/ribosomal removed, renormalised),
# cd2 (colData aligned to the columns of expr_t), tcga_group, PANEL, log_num, save_fig, fmt_p.
# ==============================================================================
suppressPackageStartupMessages({library(rstatix)})
stopifnot(identical(rownames(cd2), colnames(expr_t)) || all(cd2$barcode == colnames(expr_t)))
is_tum <- tcga_group == "Tumor"
tum <- cd2[is_tum, ]
fig_grp <- function(x) { x <- sub("^Stage ", "", x)
  ifelse(is.na(x), NA, ifelse(grepl("^IV", x), "IV", ifelse(grepl("^III", x), "III", ifelse(grepl("^II", x), "II", ifelse(grepl("^I", x), "I", NA))))) }
tum$FIGO <- factor(fig_grp(tum$figo_stage), levels = c("I", "II", "III", "IV"))
tum$Histology <- factor(dplyr::case_when(grepl("^Adenosquamous", tum$primary_diagnosis) ~ "Adenosquamous",
                                         grepl("denocarcinoma", tum$primary_diagnosis) ~ "Adenocarcinoma",
                                         grepl("quamous", tum$primary_diagnosis) ~ "Squamous", TRUE ~ NA_character_),
                        levels = c("Squamous", "Adenocarcinoma", "Adenosquamous"))
hc_call <- tum[["paper_CLIN.HPV_Hcall"]]
tum$HPV <- factor(dplyr::case_when(is.na(hc_call) ~ NA_character_, hc_call == "HPV16" ~ "HPV16", hc_call == "HPV18" ~ "HPV18",
                                   hc_call == "negative" ~ "HPV-negative", TRUE ~ "Other HPV types"),
                  levels = c("HPV16", "HPV18", "Other HPV types", "HPV-negative"))
tum$age <- as.numeric(tum$age_at_index)

# ---- Cohort descriptors (Table 3) ---------------------------------------------
tb <- function(v, nm) { t <- table(v, useNA = "ifany"); data.frame(variable = nm, level = ifelse(is.na(names(t)), "missing", names(t)), n = as.integer(t), value = NA_character_) }
desc <- rbind(
  data.frame(variable = "Primary tumours", level = "n", n = nrow(tum), value = NA_character_),
  data.frame(variable = "Age (years)", level = "mean (range)", n = NA_integer_, value = sprintf("%.1f (%d-%d)", mean(tum$age, na.rm = TRUE), min(tum$age, na.rm = TRUE), max(tum$age, na.rm = TRUE))),
  tb(tum$Histology, "Histology"), tb(tum$FIGO, "FIGO stage"), tb(tum$HPV, "HPV genotype"), tb(tum$vital_status, "Vital status"))
write.csv(desc, "results/tables/Table3_TCGA_cohort_descriptors.csv", row.names = FALSE)
print(desc)
log_num("tcga_age_mean", round(mean(tum$age, na.rm = TRUE), 1)); log_num("tcga_age_min", min(tum$age, na.rm = TRUE)); log_num("tcga_age_max", max(tum$age, na.rm = TRUE))
for (l in levels(tum$Histology)) log_num(paste0("tcga_hist_n_", l), sum(tum$Histology == l, na.rm = TRUE))
for (l in levels(tum$FIGO)) log_num(paste0("tcga_FIGO_n_", l), sum(tum$FIGO == l, na.rm = TRUE)); log_num("tcga_FIGO_missing", sum(is.na(tum$FIGO)))
for (l in levels(tum$HPV)) log_num(paste0("tcga_HPV_n_", l), sum(tum$HPV == l, na.rm = TRUE)); log_num("tcga_HPV_annotated", sum(!is.na(tum$HPV)))

# ---- Annotated heatmap: same 18-gene panel / distance / linkage / k = 2 as Fig. 5 ---
m <- expr_t[PANEL$gene, ]
z <- t(scale(t(m))); z[z > 2] <- 2; z[z < -2] <- -2
hc_cols <- hclust(dist(t(z), method = "euclidean"), method = "complete")
cl <- cutree(hc_cols, k = 2)
cl_lab <- factor(paste0("Cluster ", LETTERS[cl]))
ann_col <- data.frame(Sample = tcga_group, Cluster = cl_lab, row.names = colnames(m))
ann_col$Histology <- NA; ann_col$FIGO <- NA; ann_col$HPV <- NA
ann_col[rownames(tum), "Histology"] <- as.character(tum$Histology)
ann_col[rownames(tum), "FIGO"] <- as.character(tum$FIGO)
ann_col[rownames(tum), "HPV"] <- as.character(tum$HPV)
ann_col$Histology <- factor(ann_col$Histology, levels = levels(tum$Histology)); ann_col$FIGO <- factor(ann_col$FIGO, levels = levels(tum$FIGO)); ann_col$HPV <- factor(ann_col$HPV, levels = levels(tum$HPV))
ann_row <- data.frame(Axis = PANEL$axis, row.names = PANEL$gene)
ann_cols <- list(Sample = c(Normal = "#E69F00", Tumor = "#009E73"), Cluster = c("Cluster A" = "#CC79A7", "Cluster B" = "#56B4E9"),
                 Histology = c(Squamous = "#33a02c", Adenocarcinoma = "#6a3d9a", Adenosquamous = "#ff7f00"),
                 FIGO = c(I = "#deebf7", II = "#9ecae1", III = "#fc9272", IV = "#a50f15"),
                 HPV = c(HPV16 = "#d73027", HPV18 = "#fc8d59", `Other HPV types` = "#fee090", `HPV-negative` = "#91bfdb"))
for (ext in c("png", "tiff")) {
  f <- file.path("results/figures", paste0("Fig7_heatmap_TCGA_clinical.", ext))
  if (ext == "png") png(f, width = 3200, height = 2100, res = 300) else tiff(f, width = 3200, height = 2100, res = 300, compression = "lzw", type = "cairo")
  pheatmap::pheatmap(z, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100), clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean", clustering_method = "complete", annotation_col = ann_col, annotation_row = ann_row,
           annotation_colors = ann_cols, show_colnames = FALSE, fontsize_row = 9, treeheight_col = 30, treeheight_row = 30)
  dev.off()
}
# Association of the two clusters with each annotation (tumours only; Fisher exact, simulated p for tables larger than 2 x 2)
cl_t <- cl_lab[is_tum]
assoc <- function(a) { ok <- !is.na(a); tb <- table(cl_t[ok], droplevels(a[ok]))
  list(p = fisher.test(tb, simulate.p.value = TRUE, B = 1e5)$p.value, n = sum(ok), tb = tb) }
a_h <- assoc(tum$Histology); a_f <- assoc(tum$FIGO); a_v <- assoc(tum$HPV)
log_num("tcga_cluster_vs_histology_p", signif(a_h$p, 2)); log_num("tcga_cluster_vs_FIGO_p", signif(a_f$p, 2)); log_num("tcga_cluster_vs_HPV_p", signif(a_v$p, 2))
log_num("tcga_cluster_vs_histology_n", a_h$n); log_num("tcga_cluster_vs_FIGO_n", a_f$n); log_num("tcga_cluster_vs_HPV_n", a_v$n)
tabs <- list(histology = a_h$tb, FIGO = a_f$tb, HPV = a_v$tb)
cl_tab <- do.call(rbind, lapply(names(tabs), function(n) { d <- as.data.frame(tabs[[n]], stringsAsFactors = FALSE); names(d) <- c("cluster", "level", "n"); cbind(annotation = n, d) }))
write.csv(cl_tab, "results/tables/TableS7b_cluster_by_clinical_counts.csv", row.names = FALSE)
write.csv(data.frame(annotation = c("histology", "FIGO", "HPV"), n = c(a_h$n, a_f$n, a_v$n), fisher_p = c(a_h$p, a_f$p, a_v$p)), "results/tables/TableS7_cluster_vs_clinical_tests.csv", row.names = FALSE)
cat("\nCluster vs annotation:\n"); print(tabs); cat("p:", a_h$p, a_f$p, a_v$p, "\n")

# ---- Panel genes by histology, FIGO and HPV genotype (Kruskal-Wallis, BH across 18 genes; Dunn BH post hoc) ---
et <- expr_t[PANEL$gene, is_tum]
dc <- data.frame(t(et), check.names = FALSE); dc$Histology <- tum$Histology; dc$FIGO <- tum$FIGO; dc$HPV <- tum$HPV
long_c <- tidyr::pivot_longer(dc, all_of(PANEL$gene), names_to = "gene", values_to = "log2tpm")
kw_all <- do.call(rbind, lapply(c("Histology", "FIGO", "HPV"), function(v) {
  d <- long_c[!is.na(long_c[[v]]), ]; d$g <- droplevels(d[[v]])
  r <- d |> dplyr::group_by(gene) |> dplyr::summarise(n = dplyr::n(), KW_p = kruskal.test(log2tpm ~ g)$p.value, .groups = "drop")
  r$KW_q_BH <- p.adjust(r$KW_p, "BH"); r$factor <- v; r }))
write.csv(kw_all, "results/tables/TableS8_panel_by_clinical_KW.csv", row.names = FALSE)
for (v in c("Histology", "FIGO", "HPV")) { s <- kw_all[kw_all$factor == v, ]; log_num(paste0("tcga_panel_KW_q_lt_0.05_", v), sum(s$KW_q_BH < 0.05)); log_num(paste0("tcga_panel_KW_genes_q_lt_0.05_", v), paste(s$gene[s$KW_q_BH < 0.05], collapse = ",")) }
print(as.data.frame(kw_all |> dplyr::filter(KW_q_BH < 0.05) |> dplyr::arrange(factor, KW_q_BH)), digits = 3)
