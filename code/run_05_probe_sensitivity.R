setwd("/Users/micael/Documents/Artigo Transcriptomica Cancer Colo do Utero/Analises_v2")
sink("logs/05_probe_sensitivity.log", split = TRUE)
tryCatch({
  source("code/cells/c12_probe_rule_sensitivity.R", echo = FALSE)
  nl <- data.frame(key = names(NUM_LOG), value = vapply(NUM_LOG, function(x) paste(x, collapse = ";"), ""))
  write.csv(nl, "results/tables/numbers_logged_all.csv", row.names = FALSE); cat("\nDONE_OK\n")
}, error = function(e) { cat("\nERROR:", conditionMessage(e), "\n"); print(conditionCall(e)) })
sink()
