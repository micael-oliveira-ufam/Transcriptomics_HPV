# ==============================================================================
# Cell 06 - TCGA-CESC: GDC download (STAR - Counts workflow), sample selection,
# gene collapsing, MT-/ribosomal removal and TPM renormalisation
# ==============================================================================
suppressPackageStartupMessages({library(TCGAbiolinks); library(SummarizedExperiment)})
q <- GDCquery(project = "TCGA-CESC", data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(q, method = "api", files.per.chunk = 20, directory = "data/GDCdata")
se <- GDCprepare(q, directory = "data/GDCdata")
saveRDS(se, "data/tcga_cesc_se.rds")
cd <- as.data.frame(colData(se))
print(table(cd$sample_type))
log_num("tcga_n_samples_raw", ncol(se))
for (s in names(table(cd$sample_type))) log_num(paste0("tcga_n_", gsub(" ", "_", s)), sum(cd$sample_type == s))

# ---- Keep primary tumours (one aliquot per patient) and solid-tissue normals ----
keep <- cd$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
cd2 <- cd[keep, ]; cd2$barcode <- rownames(cd2)
cd2 <- cd2[order(cd2$sample_type, cd2$barcode), ]
cd2 <- cd2[!duplicated(paste(cd2$patient, cd2$sample_type)), ]
log_num("tcga_n_tumor_used", sum(cd2$sample_type == "Primary Tumor"))
log_num("tcga_n_normal_used", sum(cd2$sample_type == "Solid Tissue Normal"))
se2 <- se[, cd2$barcode]
tpm <- assay(se2, "tpm_unstrand")
rd  <- as.data.frame(rowData(se2))
gname <- rd$gene_name
ok <- !is.na(gname) & gname != ""
tpm <- tpm[ok, ]; gname <- gname[ok]
# one row per gene symbol: the Ensembl entry with the highest mean TPM
ord <- order(rowMeans(tpm), decreasing = TRUE)
tpm <- tpm[ord, ]; gname <- gname[ord]
dup <- duplicated(gname)
tpm <- tpm[!dup, ]; rownames(tpm) <- gname[!dup]
log_num("tcga_n_genes_collapsed", nrow(tpm))

# ---- Remove MT- and ribosomal-protein transcripts, then renormalise TPM -------
excl <- grepl(EXCLUDE_REGEX, rownames(tpm))
log_num("tcga_n_genes_excluded_MT_ribosomal", sum(excl))
log_num("tcga_n_MT_genes_present", sum(grepl("^MT-", rownames(tpm))))
log_num("tcga_median_pct_TPM_removed", round(median(colSums(tpm[excl, ]) / colSums(tpm) * 100), 2))
write.csv(data.frame(gene = rownames(tpm)[excl]), "results/tables/TableS_excluded_genes_TCGA.csv", row.names = FALSE)
tpm <- tpm[!excl, ]
tpm <- sweep(tpm, 2, colSums(tpm), "/") * 1e6       # TPM renormalised to 1e6 per sample
expr_t <- log2(tpm + 1)
log_num("tcga_n_genes_analysed", nrow(expr_t))
tcga_group <- factor(ifelse(cd2$sample_type == "Primary Tumor", "Tumor", "Normal"), levels = c("Normal", "Tumor"))
