## Global config for Oncobiome project
project_root <- "."

suppressPackageStartupMessages({
  library(here)
})

here::i_am("code/R/_config_Oncobiome.R")

dir_data_raw        <- file.path(project_root, "data", "raw")
dir_data_processed  <- file.path(project_root, "data", "processed")
dir_results         <- file.path(project_root, "results")
dir_docs            <- file.path(project_root, "docs")
dir_logs            <- file.path(project_root, "logs")

message("[config] Oncobiome config loaded. project_root = ", project_root)
