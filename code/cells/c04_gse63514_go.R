# ==============================================================================
# Cell 04 - GO Biological Process over-representation of up-regulated genes (Fig. 3)
# ==============================================================================
up_genes <- deg$gene[deg$class == "Up in cancer"]
map <- suppressMessages(bitr(rownames(expr), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
universe <- unique(map$ENTREZID)                       # all genes measured on the array
up_entrez <- unique(map$ENTREZID[map$SYMBOL %in% up_genes])
ego <- enrichGO(gene = up_entrez, universe = universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                minGSSize = 10, maxGSSize = 500, readable = TRUE)
ego_s <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)   # remove redundant terms
res <- as.data.frame(ego_s)
write.csv(as.data.frame(ego),  "results/tables/TableS3_GSE63514_GO_BP_all.csv",        row.names = FALSE)
write.csv(res,                 "results/tables/TableS3_GSE63514_GO_BP_simplified.csv", row.names = FALSE)
n_annot <- length(intersect(up_entrez, unique(AnnotationDbi::select(org.Hs.eg.db, up_entrez, "GOALL", "ENTREZID") %>%
                                              dplyr::filter(ONTOLOGYALL == "BP") %>% dplyr::pull(ENTREZID))))
log_num("go_n_input_genes", length(up_entrez)); log_num("go_n_input_with_BP", n_annot)
log_num("go_n_terms_raw", nrow(as.data.frame(ego))); log_num("go_n_terms_simplified", nrow(res))
top15 <- res %>% arrange(p.adjust) %>% slice_head(n = 15) %>%
  mutate(GeneRatioNum = Count / n_annot, Description = factor(Description, levels = rev(Description)))
write.csv(top15, "results/tables/Fig3_top15_GO_BP.csv", row.names = FALSE)
log_num("go_top15_max_padj", signif(max(top15$p.adjust), 3)); log_num("go_top1_term", as.character(top15$Description[1]))
log_num("go_top1_padj", signif(top15$p.adjust[1], 3))
p3 <- ggplot(top15, aes(GeneRatioNum, Description, size = Count, colour = p.adjust)) +
  geom_point() + scale_colour_viridis_c(option = "plasma", direction = -1, name = "Adjusted P") +
  scale_size_continuous(name = "Genes") +
  labs(x = "Gene ratio", y = "Biological processes") + theme_paper
save_fig(p3, "Fig3_GO_BP_GSE63514", 8, 5)
print(p3)
