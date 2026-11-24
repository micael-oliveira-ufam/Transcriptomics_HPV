# ==============================================================================
# Cell 01 - Environment, global parameters and helper functions
# Working directory must be the project root (the folder that contains code/).
# ==============================================================================
suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(limma)
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggrepel); library(pheatmap)
  library(RColorBrewer); library(matrixStats)
  library(clusterProfiler); library(org.Hs.eg.db)
})
set.seed(2026)                                   # fixed seed used in every script
options(timeout = 1200)                          # large GEO/GDC downloads
for (d in c("data", "results/figures", "results/tables"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Thresholds (identical to the Methods section of the manuscript) ----------
FDR_CUT <- 0.05          # Benjamini-Hochberg FDR
LFC_CUT <- 1             # |log2 fold change| threshold for the volcano / GO input

# ---- Histological groups of GSE63514 -----------------------------------------
GROUP_LEVELS <- c("Normal", "CIN1", "CIN2", "CIN3", "Cancer")
GROUP_COLS   <- c(Normal = "#4C78A8", CIN1 = "#72B7B2", CIN2 = "#F2CF5B",
                  CIN3 = "#F58518", Cancer = "#E45756")

# ---- Mitochondrially-encoded (MT-) and cytosolic ribosomal-protein transcripts
# are removed BEFORE any differential-expression / clustering step because they
# are dominated by pre-analytical (bench/tissue-quality) variation.
# MRPL*/MRPS* (nuclear-encoded mitochondrial ribosomal proteins) are retained.
EXCLUDE_REGEX <- "^MT-|^RP[SL][0-9]+[AXY]?[0-9]?L?[0-9]?$|^RPLP[0-2]$|^RPSA$"
# Pattern checked to EXCLUDE ribosomal-protein genes (RPL10A, RPS4X, RPL36AL, RPL22L1...) but KEEP kinases
# (RPS6KA1, RPS6KB1), RPL/RPS pseudogenes (RPL21P28) and mitochondrial ribosomal genes (MRPL*, MRPS*).

# ---- Registry of every number quoted in the manuscript text -------------------
NUM_LOG <- list()
log_num <- function(key, value) { NUM_LOG[[key]] <<- value; invisible(value) }

# ---- Figure writer: 300 dpi PNG + vector PDF + LZW-compressed TIFF ------------
save_fig <- function(p, name, w = 7, h = 5) {
  base <- file.path("results/figures", name)
  ggsave(paste0(base, ".png"),  p, width = w, height = h, dpi = 300, type = "cairo")
  ggsave(paste0(base, ".pdf"),  p, width = w, height = h, device = "pdf")
  ggsave(paste0(base, ".tiff"), p, width = w, height = h, dpi = 300, compression = "lzw")
  invisible(p)
}
fmt_p <- function(p) ifelse(p < 1e-3, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))
stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
theme_paper <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"))
cat("Setup OK | R", as.character(getRversion()), "| Bioconductor", as.character(BiocManager::version()), "\n")
