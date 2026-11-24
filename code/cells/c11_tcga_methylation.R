# ==============================================================================
# Cell 11 - TCGA-CESC DNA methylation (Illumina HM450, GDC SeSAMe level-3 beta values):
# tumour (01) vs solid-tissue normal (11) differentially methylated probes (DMPs), limma on M-values.
# Download is cached in data/TCGA_CESC_HM450_beta.rds (about 4 GB of raw files). Probe annotation: GEO GPL13534.
# ==============================================================================
suppressPackageStartupMessages({library(TCGAbiolinks); library(GEOquery); library(patchwork)})
rds_meth <- "data/TCGA_CESC_HM450_beta.rds"
if (!file.exists(rds_meth)) {
  filt <- list(op = "and", content = list(
    list(op = "=", content = list(field = "cases.project.project_id", value = "TCGA-CESC")),
    list(op = "=", content = list(field = "data_type", value = "Methylation Beta Value")),
    list(op = "=", content = list(field = "platform", value = "Illumina Human Methylation 450"))))
  body <- jsonlite::toJSON(list(filters = filt, size = 1000, format = "JSON", fields = "file_id,file_name,cases.samples.submitter_id,cases.samples.sample_type"), auto_unbox = TRUE)
  hits <- httr::content(httr::POST("https://api.gdc.cancer.gov/files", httr::content_type_json(), body = body), as = "parsed")$data$hits
  rs <- data.frame(file_id = vapply(hits, function(h) h$file_id, ""), file_name = vapply(hits, function(h) h$file_name, ""),
                   barcode = vapply(hits, function(h) h$cases[[1]]$samples[[1]]$submitter_id, ""), stringsAsFactors = FALSE)
  cat("GDC HM450 files:", nrow(rs), "\n"); stopifnot(!anyDuplicated(rs$barcode))
  dest <- "data/GDCdata/TCGA-CESC/DNA_Methylation/Methylation_Beta_Value"
  fl <- file.path(dest, rs$file_id, rs$file_name)
  if (!all(file.exists(fl))) {
    q_meth <- GDCquery(project = "TCGA-CESC", data.category = "DNA Methylation", data.type = "Methylation Beta Value", platform = "Illumina Human Methylation 450")
    GDCdownload(q_meth, method = "api", files.per.chunk = 10, directory = "data/GDCdata")
  }
  stopifnot(all(file.exists(fl)))
  p0 <- data.table::fread(fl[1], header = FALSE, col.names = c("probe", "beta"))$probe
  beta_all <- matrix(NA_real_, nrow = length(p0), ncol = length(fl), dimnames = list(p0, rs$barcode))
  for (i in seq_along(fl)) { d <- data.table::fread(fl[i], header = FALSE, col.names = c("probe", "beta")); beta_all[, i] <- d$beta[match(p0, d$probe)] }
  saveRDS(beta_all, rds_meth); rm(d); gc()
}
beta_all <- readRDS(rds_meth)
stype <- substr(colnames(beta_all), 14, 15)          # 01 primary tumour, 11 solid-tissue normal, 06 metastasis (not used)
log_num("meth_samples_raw", ncol(beta_all)); log_num("meth_tumour", sum(stype == "01")); log_num("meth_normal", sum(stype == "11")); log_num("meth_excluded_other_types", sum(!stype %in% c("01", "11")))
beta_all <- beta_all[, stype %in% c("01", "11")]; stype <- stype[stype %in% c("01", "11")]
if (!file.exists("data/GPL13534.soft.gz")) GEOquery::getGEO("GPL13534", destdir = "data", AnnotGPL = FALSE)
gpl <- Table(getGEO(filename = "data/GPL13534.soft.gz", AnnotGPL = FALSE)); gpl <- gpl[match(rownames(beta_all), gpl$ID), ]
s1 <- rowMeans(is.na(beta_all)) <= 0.05; b1 <- beta_all[s1, ]; g1 <- gpl[s1, ]
s2 <- !(g1$CHR %in% c("X", "Y")); b2 <- b1[s2, ]; g2 <- g1[s2, ]
snp <- (!is.na(g2$Probe_SNPs) & g2$Probe_SNPs != "") | (!is.na(g2$Probe_SNPs_10) & g2$Probe_SNPs_10 != "")
b3 <- b2[!snp, ]; g3 <- g2[!snp, ]
cnt <- c(raw = nrow(beta_all), na_le_5pct = nrow(b1), no_sex_chr = nrow(b2), no_SNP = nrow(b3)); print(cnt)
for (n in names(cnt)) log_num(paste0("meth_probes_", n), cnt[[n]])
beta2m <- function(b, off = 0.001) { b <- pmax(pmin(b, 1 - off), off); log2(b / (1 - b)) }
M <- beta2m(b3); med <- apply(M, 2, median, na.rm = TRUE); M <- sweep(M, 2, med - median(med), "-")     # median-centred M-values
grp <- factor(ifelse(stype == "01", "Tumour", "Normal"), levels = c("Normal", "Tumour"))
des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
fit <- eBayes(contrasts.fit(lmFit(M, des), makeContrasts(Tumour - Normal, levels = des)), robust = TRUE)
dmp <- topTable(fit, number = Inf, adjust.method = "BH", sort.by = "P")
dmp$status <- "Not significant"; dmp$status[dmp$adj.P.Val < 0.05 & dmp$logFC > 1] <- "Hypermethylated"; dmp$status[dmp$adj.P.Val < 0.05 & dmp$logFC < -1] <- "Hypomethylated"
dmp$probe <- rownames(dmp); gi <- match(dmp$probe, g3$ID); dmp$gene <- g3$UCSC_RefGene_Name[gi]; dmp$island_relation <- g3$Relation_to_UCSC_CpG_Island[gi]
dmp$delta_beta <- rowMeans(b3[dmp$probe, grp == "Tumour"], na.rm = TRUE) - rowMeans(b3[dmp$probe, grp == "Normal"], na.rm = TRUE)
write.csv(dmp, "results/tables/TableS12_TCGA_CESC_DMPs_all_probes.csv", row.names = FALSE)
log_num("dmp_n_tumour", sum(grp == "Tumour")); log_num("dmp_n_normal", sum(grp == "Normal"))
log_num("dmp_tested", nrow(dmp)); log_num("dmp_hyper", sum(dmp$status == "Hypermethylated")); log_num("dmp_hypo", sum(dmp$status == "Hypomethylated"))
log_num("dmp_hyper_dbeta_gt_0.2", sum(dmp$status == "Hypermethylated" & dmp$delta_beta > 0.2)); log_num("dmp_hypo_dbeta_lt_-0.2", sum(dmp$status == "Hypomethylated" & dmp$delta_beta < -0.2))
print(table(dmp$status))
# Methylation of the panel genes (promoter-region probes: TSS200/TSS1500/5'UTR/1stExon)
reg <- g3$UCSC_RefGene_Group[match(dmp$probe, g3$ID)]
prom <- grepl("TSS200|TSS1500|5'UTR|1stExon", reg)
pg <- sapply(strsplit(dmp$gene, ";"), function(x) any(x %in% PANEL$gene))
pan <- dmp[pg & prom & dmp$status != "Not significant", c("probe", "gene", "logFC", "delta_beta", "adj.P.Val", "status")]
write.csv(pan, "results/tables/TableS13_panel_gene_promoter_DMPs.csv", row.names = FALSE)
log_num("dmp_panel_promoter_sig", nrow(pan)); log_num("dmp_panel_promoter_genes", paste(sort(unique(unlist(strsplit(pan$gene, ";")))[unique(unlist(strsplit(pan$gene, ";"))) %in% PANEL$gene]), collapse = ","))
print(pan, digits = 3)

# ---- Figure 9: (A) volcano, (B) CpG-island context ------------------------------------
cols_m <- c(Hypermethylated = "#b2182b", Hypomethylated = "#2166ac", `Not significant` = "#b0b0b0"); dmp$status <- factor(dmp$status, levels = names(cols_m))
pA <- ggplot(dmp, aes(logFC, -log10(adj.P.Val), colour = status)) + geom_point(size = 0.25, alpha = 0.5, stroke = 0) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40", linewidth = 0.3) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  scale_colour_manual(values = cols_m, name = NULL, labels = function(x) paste0(sub("^Hyper.*", "Hyper", sub("^Hypo.*", "Hypo", sub("^Not.*", "NS", x))), " (n = ", format(table(dmp$status)[x], big.mark = ",", trim = TRUE), ")")) +
  labs(x = "M-value difference (tumour - normal)", y = expression(-log[10]~"adjusted p")) + guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 1)) +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
ctx <- function(x) factor(ifelse(x == "" | is.na(x), "Open sea", sub("^N_|^S_", "", x)), levels = c("Island", "Shore", "Shelf", "Open sea"))
isl <- dmp |> dplyr::mutate(ctx = ctx(island_relation)); a1 <- isl |> dplyr::filter(status != "Not significant") |> dplyr::count(status, ctx); a0 <- isl |> dplyr::count(ctx) |> dplyr::mutate(status = "All tested")
isl_all <- dplyr::bind_rows(a1 |> dplyr::mutate(status = as.character(status)), a0) |> dplyr::group_by(status) |> dplyr::mutate(pct = 100 * n / sum(n)) |> dplyr::ungroup() |>
  dplyr::mutate(status = factor(status, levels = c("All tested", "Hypermethylated", "Hypomethylated")))
write.csv(isl_all, "results/tables/TableS14_DMP_CpG_island_context.csv", row.names = FALSE)
for (s in c("All tested", "Hypermethylated", "Hypomethylated")) log_num(paste0("dmp_island_pct_", gsub(" ", "_", s)), round(isl_all$pct[isl_all$status == s & isl_all$ctx == "Island"], 1))
pB <- ggplot(isl_all, aes(status, pct, fill = ctx)) + geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(Island = "#1b7837", Shore = "#7fbf7b", Shelf = "#d9f0d3", `Open sea` = "#f7f7f7"), name = NULL) +
  labs(x = NULL, y = "Probes (%)") + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
fig9 <- pA + pB + plot_layout(widths = c(1.4, 1)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 10))
save_fig(fig9, "Fig9_methylation_TCGA_CESC", 7.2, 3.6)
