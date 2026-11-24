# Runner: executes the GSE63514 cells in order.
# NOTE: sink() interferes with GEOquery's internal connections, so logging starts after data loading.
ROOT <- getwd()  # run R from the repository root
setwd(ROOT)
run_cell <- function(f) { cat("\n########", f, "########\n"); source(file.path("code/cells", f), echo = FALSE, local = FALSE) }
run_cell("c01_setup.R")
run_cell("c02_gse63514_prepare.R")
sink("logs/01_gse63514.log", split = TRUE)
tryCatch({
  cat("c02 summary: genes collapsed =", NUM_LOG$gse63514_n_genes_collapsed, "| excluded =", NUM_LOG$gse63514_n_genes_excluded_MT_ribosomal,
      "| analysed =", NUM_LOG$gse63514_n_genes_analysed, "\n")
  for (f in c("c03_gse63514_deg.R","c04_gse63514_go.R","c05_gse63514_panel.R")) run_cell(f)
  nl <- data.frame(key = names(NUM_LOG), value = vapply(NUM_LOG, function(x) paste(x, collapse = ";"), ""))
  write.csv(nl, "results/tables/numbers_logged_GSE63514.csv", row.names = FALSE)
  cat("\nDONE_OK\n")
}, error = function(e) { cat("\nERROR:", conditionMessage(e), "\n"); print(conditionCall(e)) })
sink()
