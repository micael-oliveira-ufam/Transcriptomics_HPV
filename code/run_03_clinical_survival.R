# Runs cells 09 (clinical) and 10 (survival). Needs expr_t, cd2, tcga_group, PANEL, NUM_LOG and helpers in the session.
setwd("/Users/micael/Documents/Artigo Transcriptomica Cancer Colo do Utero/Analises_v2")
sink("logs/03_clinical_survival.log", split = TRUE)
tryCatch({
  source("code/cells/c09_tcga_clinical.R", echo = FALSE); source("code/cells/c10_tcga_survival.R", echo = FALSE)
  nl <- data.frame(key = names(NUM_LOG), value = vapply(NUM_LOG, function(x) paste(x, collapse = ";"), ""))
  write.csv(nl, "results/tables/numbers_logged_all.csv", row.names = FALSE); cat("\nDONE_OK\n")
}, error = function(e) { cat("\nERROR:", conditionMessage(e), "\n"); print(conditionCall(e)) })
sink()
