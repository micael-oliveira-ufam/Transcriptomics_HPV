# ==============================================================================
# Cell 11b - Figure 9 (methylation) with gene labels on the volcano plot.
# Re-draws Fig 9 from the tables written by cell 11 (Table S12 = all probes, S13 = panel promoter DMPs,
# S14 = CpG-island context), so it can be run without re-downloading the 4 GB of methylation data.
# Labels: (i) the panel genes with a significant promoter DMP (Table S13; bold, coloured by direction) and
# (ii) the 5 most significant hyper- and hypo-methylated probes outside the panel (grey italics).
# Requires cell 01 (save_fig, PANEL) to have been sourced in the session.
# ==============================================================================
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(dplyr)})
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)
if (!exists("PANEL")) PANEL <- data.frame(gene = c("MKI67","TOP2A","MCM2","PCNA","CDKN2A","CCND1","TP53","RB1","PIK3CA","PTEN","EGFR","ERBB2","MYC","KRT10","KRT14","HIF1A","VEGFA","CD274"))  # same 18 genes as cell 05
f12 <- "results/tables/TableS12_TCGA_CESC_DMPs_all_probes.csv"
dmp <- if (file.exists(f12)) read.csv(f12, stringsAsFactors = FALSE) else read.csv(gzfile(paste0(f12, ".gz")), stringsAsFactors = FALSE)
s13 <- read.csv("results/tables/TableS13_panel_gene_promoter_DMPs.csv", stringsAsFactors = FALSE)
isl_all <- read.csv("results/tables/TableS14_DMP_CpG_island_context.csv", stringsAsFactors = FALSE)
isl_all$status <- factor(isl_all$status, levels = c("All tested", "Hypermethylated", "Hypomethylated"))
isl_all$ctx <- factor(isl_all$ctx, levels = c("Island", "Shore", "Shelf", "Open sea"))
cols_m <- c(Hypermethylated = "#b2182b", Hypomethylated = "#2166ac", `Not significant` = "#b0b0b0")
dmp$status <- factor(dmp$status, levels = names(cols_m)); dmp$nlp <- -log10(dmp$adj.P.Val)
first_gene <- function(x) vapply(strsplit(x, ";"), function(g) { g <- unique(g[g != ""]); if (length(g)) g[1] else NA_character_ }, "")
panel_gene <- function(x) vapply(strsplit(x, ";"), function(g) { g <- unique(g[g %in% PANEL$gene]); if (length(g)) g[1] else NA_character_ }, "")
# (i) panel genes: one label per gene (probe with the smallest adjusted P among its promoter DMPs)
lab_panel <- s13 |> dplyr::mutate(label = panel_gene(gene)) |> dplyr::filter(!is.na(label)) |>
  dplyr::group_by(label) |> dplyr::slice_min(adj.P.Val, n = 1, with_ties = FALSE) |> dplyr::ungroup() |>
  dplyr::select(probe, label) |> dplyr::left_join(dmp[, c("probe", "logFC", "nlp", "status")], by = "probe") |> dplyr::mutate(grp = "Panel gene")
# (ii) top 5 hyper + top 5 hypo (by adjusted P) among non-panel genes
lab_top <- dmp |> dplyr::filter(status != "Not significant", is.na(panel_gene(gene)), !is.na(first_gene(gene))) |>
  dplyr::mutate(label = first_gene(gene)) |> dplyr::group_by(status, label) |> dplyr::slice_min(adj.P.Val, n = 1, with_ties = FALSE) |> dplyr::ungroup() |>  # one probe per gene
  dplyr::group_by(status) |> dplyr::slice_min(adj.P.Val, n = 5, with_ties = FALSE) |> dplyr::ungroup() |>
  dplyr::mutate(grp = "Top probe") |> dplyr::select(probe, label, logFC, nlp, status, grp)
lab <- dplyr::bind_rows(lab_panel, lab_top)
write.csv(lab, "results/tables/Fig9A_labelled_genes.csv", row.names = FALSE)
pA <- ggplot(dmp, aes(logFC, nlp, colour = status)) + geom_point(size = 0.25, alpha = 0.5, stroke = 0) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_point(data = lab, shape = 21, fill = NA, size = 1.6, stroke = 0.4, colour = "black", show.legend = FALSE) +
  geom_text_repel(data = lab |> dplyr::filter(grp == "Top probe"), aes(label = label), size = 2, fontface = "italic", colour = "grey25",
                  segment.size = 0.2, segment.colour = "grey50", min.segment.length = 0, max.overlaps = Inf, box.padding = 0.35, seed = 2026, bg.color = "white", bg.r = 0.12, nudge_y = 1.2, show.legend = FALSE) +
  geom_text_repel(data = lab |> dplyr::filter(grp == "Panel gene"), aes(label = label), size = 2.4, fontface = "bold.italic",
                  segment.size = 0.25, segment.colour = "black", min.segment.length = 0, max.overlaps = Inf, box.padding = 0.6, point.padding = 0.3, seed = 2026, bg.color = "white", bg.r = 0.15, nudge_x = ifelse(lab$logFC[lab$grp == "Panel gene"] > 0, 1.6, -1.6), nudge_y = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = cols_m, name = NULL,
                      labels = function(x) paste0(sub("^Hyper.*", "Hyper", sub("^Hypo.*", "Hypo", sub("^Not.*", "NS", x))), " (n = ", format(table(dmp$status)[x], big.mark = ",", trim = TRUE), ")")) +
  labs(x = "M-value difference (tumour - normal)", y = expression(-log[10]~"adjusted p")) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) + theme_classic(base_size = 9) + theme(legend.position = "bottom")
pB <- ggplot(isl_all, aes(status, pct, fill = ctx)) + geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(Island = "#1b7837", Shore = "#7fbf7b", Shelf = "#d9f0d3", `Open sea` = "#f7f7f7"), name = NULL) +
  labs(x = NULL, y = "Probes (%)") + theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")
fig9 <- pA + pB + plot_layout(widths = c(1.4, 1)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 10))
save_fig(fig9, "Fig9_methylation_TCGA_CESC", 7.2, 3.6)
print(as.data.frame(lab[, c("label", "grp", "status", "logFC")]), digits = 3)
