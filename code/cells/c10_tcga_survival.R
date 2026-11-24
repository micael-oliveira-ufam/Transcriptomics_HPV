# ==============================================================================
# Cell 10 - TCGA-CESC overall survival: CD274 (pre-specified) and the 18-gene panel + checkpoint/TME genes.
# Event = death; time = days to death, otherwise days to last follow-up (months = days / 30.4375).
# Cox models use z-scored log2 TPM (HR per +1 SD); BH within each model family; Kaplan-Meier on the median split.
# Needs c09 objects: tum (tumour clinical table with FIGO/Histology/age), expr_t, tcga_group, PANEL.
# ==============================================================================
suppressPackageStartupMessages({library(survival); library(survminer); library(patchwork)})
is_tum <- tcga_group == "Tumor"
sv <- data.frame(time_d = ifelse(tum$vital_status == "Dead", as.numeric(tum$days_to_death), as.numeric(tum$days_to_last_follow_up)),
                 event = as.integer(tum$vital_status == "Dead"), age = tum$age,
                 stage = factor(ifelse(tum$FIGO %in% c("I", "II"), "I-II", ifelse(tum$FIGO %in% c("III", "IV"), "III-IV", NA)), levels = c("I-II", "III-IV")),
                 hist = factor(ifelse(tum$Histology == "Squamous", "SCC", ifelse(is.na(tum$Histology), NA, "Non-SCC")), levels = c("SCC", "Non-SCC")),
                 row.names = rownames(tum))
sv$time_m <- sv$time_d / 30.4375
ok <- !is.na(sv$time_d) & sv$time_d > 0
log_num("surv_n_tumours", nrow(sv)); log_num("surv_n_used", sum(ok)); log_num("surv_events", sum(sv$event[ok])); log_num("surv_n_excluded", sum(!ok))
log_num("surv_median_followup_months", round(median(sv$time_m[ok & sv$event == 0]), 1))

genes_s <- unique(c(PANEL$gene, "CD274", "PDCD1", "CTLA4", "FOXP3", "CD8A"))
stopifnot(all(genes_s %in% rownames(expr_t)))
Xs <- t(scale(t(expr_t[genes_s, is_tum])))               # z-score across tumours
dsv <- cbind(sv[, c("time_m", "event", "age", "stage", "hist")], t(Xs))[ok, ]
uni <- do.call(rbind, lapply(genes_s, function(g) {
  s <- summary(coxph(as.formula(paste0("Surv(time_m, event) ~ `", g, "`")), data = dsv))
  data.frame(gene = g, HR = s$conf.int[1, "exp(coef)"], lo = s$conf.int[1, "lower .95"], hi = s$conf.int[1, "upper .95"], p = s$coefficients[1, "Pr(>|z|)"]) }))
uni$q_BH <- p.adjust(uni$p, "BH"); uni <- uni[order(uni$p), ]
stopifnot(all(uni$lo <= uni$HR & uni$HR <= uni$hi))
write.csv(uni, "results/tables/TableS9_Cox_univariable_HR_per_SD.csv", row.names = FALSE)
log_num("cox_genes_tested", nrow(uni)); log_num("cox_genes_q_lt_0.05", sum(uni$q_BH < 0.05)); log_num("cox_genes_p_lt_0.05", sum(uni$p < 0.05))
for (i in 1:5) log_num(paste0("cox_top", i), paste0(uni$gene[i], " HR ", round(uni$HR[i], 2), " (", round(uni$lo[i], 2), "-", round(uni$hi[i], 2), "); p = ", signif(uni$p[i], 2), "; q = ", signif(uni$q_BH[i], 2)))
log_num("cox_CD274_uni", paste0("HR ", round(uni$HR[uni$gene == "CD274"], 2), " (", round(uni$lo[uni$gene == "CD274"], 2), "-", round(uni$hi[uni$gene == "CD274"], 2), "); p = ", signif(uni$p[uni$gene == "CD274"], 2)))
print(uni, digits = 3)

# Adjusted model for CD274: age, FIGO group (I-II vs III-IV) and histology (SCC vs non-SCC); complete cases
dadj <- dsv[complete.cases(dsv[, c("age", "stage", "hist", "CD274")]), ]
fit_adj <- coxph(Surv(time_m, event) ~ CD274 + age + stage + hist, data = dadj)
sa <- summary(fit_adj)
adj <- data.frame(term = rownames(sa$coefficients), HR = sa$conf.int[, "exp(coef)"], lo = sa$conf.int[, "lower .95"], hi = sa$conf.int[, "upper .95"], p = sa$coefficients[, "Pr(>|z|)"])
write.csv(adj, "results/tables/TableS10_Cox_CD274_adjusted.csv", row.names = FALSE)
log_num("cox_adj_n", nrow(dadj)); log_num("cox_adj_events", sum(dadj$event))
log_num("cox_adj_CD274", paste0("HR ", round(adj$HR[1], 2), " (", round(adj$lo[1], 2), "-", round(adj$hi[1], 2), "); p = ", signif(adj$p[1], 2)))
log_num("cox_adj_stage", paste0("HR ", round(adj$HR[3], 2), " (", round(adj$lo[3], 2), "-", round(adj$hi[3], 2), "); p = ", signif(adj$p[3], 2)))
ph <- cox.zph(fit_adj); log_num("cox_adj_PH_global_p", signif(ph$table["GLOBAL", "p"], 2))
print(adj, digits = 3); print(ph)

# Kaplan-Meier, median split, log-rank (all genes, BH)
km_p <- sapply(genes_s, function(g) { grp <- ifelse(dsv[[g]] >= median(dsv[[g]]), "High", "Low"); 1 - pchisq(survdiff(Surv(dsv$time_m, dsv$event) ~ grp)$chisq, 1) })
km_tab <- data.frame(gene = genes_s, logrank_p = km_p, q_BH = p.adjust(km_p, "BH")); write.csv(km_tab, "results/tables/TableS11_KM_logrank.csv", row.names = FALSE)
log_num("km_CD274_logrank_p", signif(km_p["CD274"], 2)); log_num("km_genes_q_lt_0.05", sum(km_tab$q_BH < 0.05))

# ---- Figure 6: (A) KM CD274, (B) forest of univariable HR per SD -----------------
grp_cd <- factor(ifelse(dsv$CD274 >= median(dsv$CD274), "High", "Low"), levels = c("High", "Low"))
dk <- data.frame(dsv, grp_cd); fit_cd <- survfit(Surv(time_m, event) ~ grp_cd, data = dk)
pA <- ggsurvplot(fit_cd, data = dk, palette = c("#d73027", "#4575b4"), conf.int = FALSE, censor.size = 2, xlab = "Time (months)", ylab = "Overall survival",
                 legend.title = "", legend.labs = c("CD274 high", "CD274 low"), legend = c(0.75, 0.9), ggtheme = theme_classic(base_size = 9), size = 0.6)$plot +
  annotate("text", x = 0, y = 0.08, hjust = 0, size = 2.6, label = paste0("log-rank p = ", fmt_p(km_p["CD274"]), "\nn = ", nrow(dsv), ", deaths = ", sum(dsv$event)))
fp <- uni; fp$gene <- factor(fp$gene, levels = rev(fp$gene)); fp$sig <- ifelse(fp$q_BH < 0.05, "q < 0.05", "q >= 0.05")
pB <- ggplot(fp, aes(HR, gene, colour = sig)) + geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.25, linewidth = 0.4) + scale_x_log10() +
  scale_colour_manual(values = c("q < 0.05" = "#d73027", "q >= 0.05" = "#555555"), name = NULL) +
  labs(x = "Hazard ratio per +1 SD (95% CI)", y = NULL) + theme_classic(base_size = 9) + theme(axis.text.y = element_text(face = "italic", size = 6), legend.position = "bottom")
fig8 <- pA + pB + plot_layout(widths = c(1, 1.1)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 10))
save_fig(fig8, "Fig8_survival_TCGA_CESC", 7.2, 4.6)
