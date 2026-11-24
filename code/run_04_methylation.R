# Runs cell 11 (methylation). Downloads about 4 GB from GDC on the first run (cached afterwards).
setwd("/Users/micael/Documents/Artigo Transcriptomica Cancer Colo do Utero/Analises_v2")
sink("logs/04_methylation.log", split = TRUE)
tryCatch({
  source("code/cells/c11_tcga_methylation.R", echo = FALSE)
  nl <- data.frame(key = names(NUM_LOG), value = vapply(NUM_LOG, function(x) paste(x, collapse = ";"), ""))
  write.csv(nl, "results/tables/numbers_logged_all.csv", row.names = FALSE); cat("\nDONE_OK\n")
}, error = function(e) { cat("\nERROR:", conditionMessage(e), "\n"); print(conditionCall(e)) })
sink()
