# Runs only the TCGA volcano cell (needs expr_t, tcga_group, NUM_LOG; deg from run_01 if available, for the concordance with GSE63514).
ROOT <- getwd()  # run R from the repository root
sink("logs/02c_tcga_volcano.log", split = TRUE)
tryCatch({
  source("code/cells/c08_tcga_volcano.R", echo = FALSE)
  nl <- data.frame(key = names(NUM_LOG), value = vapply(NUM_LOG, function(x) paste(x, collapse = ";"), ""))
  write.csv(nl, "results/tables/numbers_logged_all.csv", row.names = FALSE); cat("\nDONE_OK\n")
}, error = function(e) { cat("\nERROR:", conditionMessage(e), "\n"); print(conditionCall(e)) })
sink()
