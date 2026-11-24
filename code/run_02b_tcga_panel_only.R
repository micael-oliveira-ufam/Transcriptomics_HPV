# Re-runs only the TCGA panel/heatmap cell (needs expr_t, tcga_group, PANEL, NUM_LOG in the R session).
setwd("/Users/micael/Documents/Artigo Transcriptomica Cancer Colo do Utero/Analises_v2")
sink("logs/02b_tcga_panel.log", split = TRUE)
tryCatch({
  source("code/cells/c07_tcga_panel.R", echo = FALSE)
  nl <- data.frame(key = names(NUM_LOG), value = vapply(NUM_LOG, function(x) paste(x, collapse = ";"), ""))
  write.csv(nl, "results/tables/numbers_logged_all.csv", row.names = FALSE); cat("\nDONE_OK\n")
}, error = function(e) { cat("\nERROR:", conditionMessage(e), "\n"); print(conditionCall(e)) })
sink()
