# Runner for TCGA-CESC. Logging starts after the GDC download (sink interferes with downloads).
ROOT <- getwd()  # run R from the repository root
setwd(ROOT)
run_cell <- function(f) { cat("\n########", f, "########\n"); source(file.path("code/cells", f), echo = FALSE, local = FALSE) }
# PANEL is defined in c05; recreate numbers registry state from the GSE63514 session if present
if (!exists("PANEL")) { run_cell("c01_setup.R") }
if (!exists("PANEL")) stop("Run run_01_gse63514.R first in this R session (PANEL object needed).")
run_cell("c06_tcga_prepare.R")
sink("logs/02_tcga.log", split = TRUE)
tryCatch({
  cat("tcga summary:", NUM_LOG$tcga_n_tumor_used, "tumours,", NUM_LOG$tcga_n_normal_used, "normals,", NUM_LOG$tcga_n_genes_analysed, "genes\n")
  run_cell("c07_tcga_panel.R")
  nl <- data.frame(key = names(NUM_LOG), value = vapply(NUM_LOG, function(x) paste(x, collapse = ";"), ""))
  write.csv(nl, "results/tables/numbers_logged_all.csv", row.names = FALSE)
  cat("\nDONE_OK\n")
}, error = function(e) { cat("\nERROR:", conditionMessage(e), "\n"); print(conditionCall(e)) })
sink()
