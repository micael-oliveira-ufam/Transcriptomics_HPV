# Runs the complete analysis from the repository root: R -f code/run_all.R   (or source() it in RStudio)
# Order matters: run_01 defines PANEL and NUM_LOG, which run_02 needs in the same R session.
stopifnot(file.exists("code/cells/c01_setup.R"))
dir.create("logs", showWarnings = FALSE)
source("code/run_01_gse63514.R")
source("code/run_02_tcga.R")
