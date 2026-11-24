# ==============================================================================
# Cell 12 - Sensitivity of the 18-transcript panel (GSE63514) to the probe-collapsing rule
# Main analysis: highest mean probe per gene. Alternative: highest-variance probe per gene (rule used in the earlier draft for another series).
# Needs from the session: eset, expr_probe, keep, sym, group, GROUP_LEVELS, PANEL, EXCLUDE_REGEX, stats (from cell 05), log_num.
# ==============================================================================
cat("Sample-characteristic fields available in the GSE63514 series matrix:\n"); print(grep("ch1$", colnames(pData(eset)), value = TRUE))
log_num("gse63514_characteristic_fields", paste(sub(":ch1$", "", grep("ch1$", colnames(pData(eset)), value = TRUE)), collapse = "; "))
# Objects below are rebuilt here because later cells (c10/c11) may have overwritten names such as keep/group in the live session
sym  <- trimws(fData(eset)[["Gene symbol"]]); keep <- !is.na(sym) & sym != "" & !grepl("///", sym, fixed = TRUE)
tt <- pData(eset)[["tissue type:ch1"]]
group <- factor(dplyr::case_when(grepl("normal", tt) ~ "Normal", grepl("low grade", tt) ~ "CIN1", grepl("moderate grade", tt) ~ "CIN2",
                                 grepl("high grade", tt) ~ "CIN3", grepl("cancer", tt) ~ "Cancer"), levels = GROUP_LEVELS)
stopifnot(length(keep) == nrow(expr_probe), length(group) == ncol(expr_probe), all(c("gene","p_anova","tukey_Cancer_vs_Normal") %in% names(stats)))
pt <- tibble(probe = rownames(expr_probe)[keep], symbol = sym[keep], v = matrixStats::rowVars(expr_probe[keep, ]))
bv <- pt %>% group_by(symbol) %>% slice_max(v, n = 1, with_ties = FALSE) %>% ungroup()
ev <- expr_probe[bv$probe, ]; rownames(ev) <- bv$symbol; ev <- ev[!grepl(EXCLUDE_REGEX, rownames(ev)), ]
bm <- tibble(probe = rownames(expr_probe)[keep], symbol = sym[keep], m = rowMeans(expr_probe[keep, ])) %>% group_by(symbol) %>% slice_max(m, n = 1, with_ties = FALSE) %>% ungroup()
same_probe <- bv$probe %in% bm$probe[match(bv$symbol, bm$symbol)]
log_num("probe_rule_genes_same_probe_pct", round(100 * mean(bv$probe == bm$probe[match(bv$symbol, bm$symbol)]), 1))
sens <- lapply(PANEL$gene, function(g) {
  y <- ev[g, ]; a <- aov(y ~ group); tk <- TukeyHSD(a)$group
  tibble(gene = g, probe_highest_mean = bm$probe[bm$symbol == g], probe_highest_variance = bv$probe[bv$symbol == g],
         anova_p_mean = stats$p_anova[stats$gene == g], anova_p_variance = summary(a)[[1]][["Pr(>F)"]][1],
         tukey_Cancer_vs_Normal_mean = stats$tukey_Cancer_vs_Normal[stats$gene == g], tukey_Cancer_vs_Normal_variance = tk["Cancer-Normal", "p adj"]) }) %>% bind_rows()
sens$bonf_sig_mean <- sens$anova_p_mean < 0.05 / 18; sens$bonf_sig_variance <- sens$anova_p_variance < 0.05 / 18
write.csv(sens, "results/tables/TableS15_panel_probe_rule_sensitivity.csv", row.names = FALSE)
log_num("probe_rule_panel_same_probe", sum(sens$probe_highest_mean == sens$probe_highest_variance))
log_num("probe_rule_bonferroni_sig_mean", sum(sens$bonf_sig_mean)); log_num("probe_rule_bonferroni_sig_variance", sum(sens$bonf_sig_variance))
log_num("probe_rule_bonferroni_changed", paste(sens$gene[sens$bonf_sig_mean != sens$bonf_sig_variance], collapse = ","))
print(as.data.frame(sens %>% mutate(across(where(is.numeric), ~ signif(., 3)))))
