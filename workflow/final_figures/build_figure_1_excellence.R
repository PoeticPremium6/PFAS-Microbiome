suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})

args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop("Repository root is required.", call.=FALSE)

REV <- normalizePath(args[[1]], mustWork=TRUE)
OUT <- file.path(REV, "submission", "figure_builds", "Figure_1")
FINAL <- file.path(OUT, "figures")
PREVIEW <- file.path(OUT, "previews")
LOGS <- file.path(OUT, "reports")
FIGSRC <- file.path(OUT, "source_data")
PANELDIR <- file.path(OUT, "panels")
SCRIPT_DIR <- file.path(REV, "scripts", "final_figures", "figure1_assets")

for (d in c(OUT, FINAL, PREVIEW, LOGS, FIGSRC, PANELDIR)) {
  dir.create(d, recursive=TRUE, showWarnings=FALSE)
}

Sys.setenv(
  PFAS_REBUILD_MAIN = FINAL,
  PFAS_REBUILD_PREVIEWS = PREVIEW,
  PFAS_REBUILD_LOGS = LOGS,
  PFAS_REBUILD_SOURCE_DATA = FIGSRC,
  PFAS_REBUILD_PANELS = PANELDIR,
  PFAS_REBUILD_SCRIPT_DIR = SCRIPT_DIR
)

source(file.path(SCRIPT_DIR, "phase14_figure1_v13_excellence_closeout_v2.R"), local=TRUE)
