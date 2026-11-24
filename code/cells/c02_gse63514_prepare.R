# ==============================================================================
# Cell 02 - GSE63514: download, group parsing, probe-to-gene collapsing
# Affymetrix Human Genome U133 Plus 2.0 (GPL570); den Boon et al., PNAS 2015
# ==============================================================================
# Cached ExpressionSet is used when present (offline-safe); otherwise it is downloaded from GEO and cached.
if (file.exists("data/gse63514_eset.rds")) {
  eset <- readRDS("data/gse63514_eset.rds")
} else {
  eset <- getGEO("GSE63514", GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data")[[1]]
  saveRDS(eset, "data/gse63514_eset.rds")
}
expr_probe <- exprs(eset)
log_num("gse63514_n_probes_raw", nrow(expr_probe))
log_num("gse63514_max_intensity", signif(max(expr_probe, na.rm = TRUE), 5))

# Scale check (Methods 2.3): log2(x + 1) only if the matrix is not already on log scale
if (max(expr_probe, na.rm = TRUE) > 100) expr_probe <- log2(expr_probe + 1)

# ---- Histological group from the sample characteristics ----------------------
pd <- pData(eset)
tt <- pd[["tissue type:ch1"]]
group <- factor(dplyr::case_when(
  grepl("normal",         tt) ~ "Normal",
  grepl("low grade",      tt) ~ "CIN1",
  grepl("moderate grade", tt) ~ "CIN2",
  grepl("high grade",     tt) ~ "CIN3",
  grepl("cancer",         tt) ~ "Cancer"), levels = GROUP_LEVELS)
stopifnot(!anyNA(group))
print(table(group))
for (g in GROUP_LEVELS) log_num(paste0("gse63514_n_", g), sum(group == g))

# Dissection method (technical covariate used only in a sensitivity analysis)
dissection <- factor(ifelse(grepl("whole", pd[["dissection:ch1"]]), "whole_section", "laser_capture"))
print(table(group, dissection))

# ---- Probe -> HGNC symbol; multi-gene ("///") and unannotated probes dropped ---
sym  <- trimws(fData(eset)[["Gene symbol"]])
keep <- !is.na(sym) & sym != "" & !grepl("///", sym, fixed = TRUE)
log_num("gse63514_n_probes_dropped_no_or_multi_symbol", sum(!keep))
probe_tab <- tibble(probe = rownames(expr_probe)[keep], symbol = sym[keep],
                    mean_expr = rowMeans(expr_probe[keep, ]))
# One probe per gene: the probe with the highest mean expression across samples
best <- probe_tab %>% group_by(symbol) %>% slice_max(mean_expr, n = 1, with_ties = FALSE) %>% ungroup()
expr_all <- expr_probe[best$probe, ]; rownames(expr_all) <- best$symbol
log_num("gse63514_n_genes_collapsed", nrow(expr_all))

# ---- Remove MT- and ribosomal-protein transcripts before DEG / clustering ------
excl <- grepl(EXCLUDE_REGEX, rownames(expr_all))
write.csv(data.frame(gene = rownames(expr_all)[excl]),
          "results/tables/TableS_excluded_genes_GSE63514.csv", row.names = FALSE)
log_num("gse63514_n_genes_excluded_MT_ribosomal", sum(excl))
log_num("gse63514_n_MT_genes_present", sum(grepl("^MT-", rownames(expr_all))))
expr <- expr_all[!excl, ]
log_num("gse63514_n_genes_analysed", nrow(expr))
cat(sprintf("GSE63514: %d genes after collapsing; %d MT-/ribosomal removed; %d analysed\n",
            nrow(expr_all), sum(excl), nrow(expr)))
