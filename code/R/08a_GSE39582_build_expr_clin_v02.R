## 08a_GSE39582_build_expr_clin_v02.R
## Build proper expression + clinical objects for GSE39582 (v02)

source("./code/R/_config_Oncobiome.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c("GEOquery", "Biobase", "dplyr")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(GEOquery)
library(Biobase)
library(dplyr)

log_file <- file.path(
  dir_logs,
  paste0("08a_GSE39582_build_expr_clin_v02_",
         format(Sys.time(), "%Y%m%d_%H%M%S"),
         ".log")
)

log_message <- function(...) {
  msg <- paste0("[", Sys.time(), "] ", paste0(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_message("08a_GSE39582_build_expr_clin_v02.R started.")

external_proc_dir <- file.path(dir_data_processed, "external", "GSE39582")
dir.create(external_proc_dir, recursive = TRUE, showWarnings = FALSE)

expr_rds_v02 <- file.path(external_proc_dir, "GSE39582_expr_symbol_log2_v02.rds")
clin_rds_v02 <- file.path(external_proc_dir, "GSE39582_clinical_v02.rds")

# -----------------------------------------------------------
# Load GSE39582 from GEO
# -----------------------------------------------------------

log_message("Downloading/loading GSE39582 via GEOquery...")
gse_list <- getGEO("GSE39582", GSEMatrix = TRUE)
if (is.list(gse_list)) {
  sizes <- sapply(gse_list, ncol)
  eset <- gse_list[[which.max(sizes)]]
} else {
  eset <- gse_list
}

expr_raw <- exprs(eset)  # probes x samples
pheno    <- pData(eset)
fd       <- fData(eset)

log_message(paste0("expr_raw dim: ", nrow(expr_raw), " probes x ",
                   ncol(expr_raw), " samples."))

# -----------------------------------------------------------
# Map probes -> gene symbols using featureData
# -----------------------------------------------------------

symbol_cols <- grep("symbol", colnames(fd), ignore.case = TRUE, value = TRUE)
if (length(symbol_cols) == 0) {
  stop("No column containing 'symbol' found in featureData.")
}
sym_col <- symbol_cols[1]
log_message("Using featureData column for SYMBOL: ", sym_col)

gene_symbols <- as.character(fd[rownames(expr_raw), sym_col])

valid <- !is.na(gene_symbols) & gene_symbols != ""
expr_raw_valid     <- expr_raw[valid, , drop = FALSE]
gene_symbols_valid <- gene_symbols[valid]

log_message(paste0("Mapped ", sum(valid), " probes to non-NA gene symbols."))

expr_sum <- rowsum(expr_raw_valid, group = gene_symbols_valid)
counts   <- table(gene_symbols_valid)
expr_symbol <- sweep(expr_sum, 1, as.numeric(counts[rownames(expr_sum)]), "/")

max_val <- max(expr_symbol, na.rm = TRUE)
if (max_val > 100) {
  log_message("Max expression > 100; applying log2(x+1).")
  expr_symbol <- log2(expr_symbol + 1)
} else {
  log_message(paste0("Expression appears already log-scale (max ≈ ",
                     round(max_val, 2), ")."))
}

log_message(paste0("expr_symbol dim: ",
                   nrow(expr_symbol), " genes x ",
                   ncol(expr_symbol), " samples."))
log_message("Example colnames(expr_symbol): ",
            paste(head(colnames(expr_symbol), 5), collapse = ", "))

saveRDS(expr_symbol, expr_rds_v02)
log_message("Saved expression to: ", expr_rds_v02)

# -----------------------------------------------------------
# Clinical: OS_time / OS_event / stage_simplified
# -----------------------------------------------------------

time_col  <- "os.delay (months):ch1"
event_col <- "os.event:ch1"
stage_col <- "tnm.stage:ch1"

if (!all(c(time_col, event_col, stage_col) %in% colnames(pheno))) {
  stop("Expected clinical columns not all found in pheno.")
}

OS_time  <- suppressWarnings(as.numeric(as.character(pheno[[time_col]])))
OS_event <- suppressWarnings(as.numeric(as.character(pheno[[event_col]])))
stage_raw_num <- suppressWarnings(as.numeric(as.character(pheno[[stage_col]])))

stage_simplified <- dplyr::case_when(
  stage_raw_num %in% c(1, 2) ~ "I/II",
  stage_raw_num %in% c(3, 4) ~ "III/IV",
  TRUE ~ NA_character_
)

clin2 <- data.frame(
  sample = rownames(pheno),
  OS_time = OS_time,
  OS_event = OS_event,
  stage_simplified = stage_simplified,
  stringsAsFactors = FALSE
)
rownames(clin2) <- clin2$sample

log_message(paste0("Clinical rows: ", nrow(clin2)))
log_message("table(stage_simplified):")
capture.output(
  print(table(clin2$stage_simplified, useNA = "ifany")),
  file = log_file, append = TRUE
)

saveRDS(clin2, clin_rds_v02)
log_message("Saved clinical to: ", clin_rds_v02)

log_message("08a_GSE39582_build_expr_clin_v02.R finished.")
