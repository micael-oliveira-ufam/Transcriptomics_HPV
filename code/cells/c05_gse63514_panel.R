# ==============================================================================
# Cell 05 - 18-gene panel along cervical progression (Fig. 4)
# One-way ANOVA + Tukey HSD post hoc + linear trend test over ordered grades
# ==============================================================================
PANEL <- tibble::tribble(
  ~gene, ~axis,
  "MKI67","Proliferation / DNA replication", "TOP2A","Proliferation / DNA replication",
  "MCM2","Proliferation / DNA replication",  "PCNA","Proliferation / DNA replication",
  "CDKN2A","HPV / cell cycle", "CCND1","HPV / cell cycle", "TP53","HPV / cell cycle", "RB1","HPV / cell cycle",
  "PIK3CA","Oncogenic signalling", "PTEN","Oncogenic signalling", "EGFR","Oncogenic signalling",
  "ERBB2","Oncogenic signalling", "MYC","Oncogenic signalling",
  "KRT10","Squamous differentiation", "KRT14","Squamous differentiation",
  "HIF1A","Hypoxia / angiogenesis / immune evasion", "VEGFA","Hypoxia / angiogenesis / immune evasion",
  "CD274","Hypoxia / angiogenesis / immune evasion")
stopifnot(all(PANEL$gene %in% rownames(expr)))
score <- as.numeric(group)                        # ordinal score 1..5

stats <- lapply(PANEL$gene, function(g) {
  y  <- expr[g, ]
  a  <- aov(y ~ group); p_anova <- summary(a)[[1]][["Pr(>F)"]][1]
  tk <- TukeyHSD(a)$group
  tr <- summary(lm(y ~ score))$coefficients["score", ]
  med <- tapply(y, group, median)
  tibble(gene = g, p_anova = p_anova,
         tukey_CIN1_vs_Normal = tk["CIN1-Normal", "p adj"], tukey_CIN2_vs_Normal = tk["CIN2-Normal", "p adj"],
         tukey_CIN3_vs_Normal = tk["CIN3-Normal", "p adj"], tukey_Cancer_vs_Normal = tk["Cancer-Normal", "p adj"],
         tukey_Cancer_vs_CIN3 = tk["Cancer-CIN3", "p adj"],
         trend_slope_per_grade = tr["Estimate"], trend_p = tr["Pr(>|t|)"],
         spearman_rho = suppressWarnings(cor(y, score, method = "spearman")),
         med_Normal = med["Normal"], med_CIN1 = med["CIN1"], med_CIN2 = med["CIN2"], med_CIN3 = med["CIN3"], med_Cancer = med["Cancer"])
}) %>% bind_rows() %>% mutate(fdr_anova = p.adjust(p_anova, "BH"), fdr_trend = p.adjust(trend_p, "BH")) %>%
  left_join(PANEL, by = "gene") %>% relocate(axis, .after = gene)
write.csv(stats, "results/tables/TableS4_panel_ANOVA_Tukey_trend.csv", row.names = FALSE)
for (i in seq_len(nrow(stats))) {
  log_num(paste0("panel_", stats$gene[i], "_anova_p"), signif(stats$p_anova[i], 3))
  log_num(paste0("panel_", stats$gene[i], "_trend_slope"), round(stats$trend_slope_per_grade[i], 3))
}
log_num("panel_bonferroni_threshold", round(0.05 / nrow(PANEL), 4))
print(as.data.frame(stats %>% mutate(across(where(is.numeric), ~ signif(., 3))) %>%
                      dplyr::select(gene, p_anova, fdr_anova, starts_with("tukey"), trend_slope_per_grade, trend_p)))

# ---- Plot data, per-panel annotation (ANOVA P + stars vs Normal by Tukey HSD) ----
long <- as.data.frame(t(expr[PANEL$gene, ])) %>% mutate(group = group) %>%
  pivot_longer(-group, names_to = "gene", values_to = "expr") %>%
  mutate(gene = factor(gene, levels = PANEL$gene))
ann <- stats %>% dplyr::select(gene, p_anova, starts_with("tukey_") ) %>%
  pivot_longer(c(tukey_CIN1_vs_Normal, tukey_CIN2_vs_Normal, tukey_CIN3_vs_Normal, tukey_Cancer_vs_Normal),
               names_to = "cmp", values_to = "p") %>%
  mutate(group = factor(sub("tukey_(.*)_vs_Normal", "\\1", cmp), levels = GROUP_LEVELS), lab = stars(p)) %>%
  left_join(long %>% group_by(gene) %>% summarise(ymax = max(expr), yr = diff(range(expr))), by = "gene") %>%
  mutate(y = ymax + 0.04 * yr, gene = factor(gene, levels = PANEL$gene))
hdr <- stats %>% mutate(gene = factor(gene, levels = PANEL$gene),
                        lab = sprintf("ANOVA P = %s; Cancer vs CIN3: %s", fmt_p(p_anova), stars(tukey_Cancer_vs_CIN3))) %>%
  left_join(long %>% group_by(gene) %>% summarise(ymax = max(expr), yr = diff(range(expr))), by = "gene") %>%
  mutate(y = ymax + 0.2 * yr)
p4 <- ggplot(long, aes(group, expr, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.5) +
  geom_text(data = ann, aes(x = group, y = y, label = lab), inherit.aes = FALSE, size = 3) +
  geom_text(data = hdr, aes(x = 3, y = y, label = lab), inherit.aes = FALSE, size = 2.5) +
  facet_wrap(~ gene, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = GROUP_COLS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = NULL, y = expression(log[2]~"expression")) + theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p4, "Fig4_panel_GSE63514", 13, 7.5)
print(p4)
