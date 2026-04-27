## 03_TCGA_download_preprocess.R
## Download / load TCGA COAD+READ, normalize to gene symbols, build clinical table

source("./code/R/_config_Oncobiome.R")

cran_pkgs <- c("dplyr")
bioc_pkgs <- c("TCGAbiolinks", "SummarizedExperiment", "edgeR", "limma")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(dplyr)

log_file <- file.path(dir_logs, paste0("03_TCGA_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_message <- function(...) {
  msg <- paste0("[", Sys.time(), "] ", paste0(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_message("03_TCGA_download_preprocess.R started.")

se_rds      <- file.path(dir_data_processed, "tcga_expression", "tcga_crc_SE_v01.rds")
expr_rds    <- file.path(dir_data_processed, "tcga_expression", "tcga_crc_expr_v01.rds")
clin_rds    <- file.path(dir_data_processed, "tcga_clinical",  "tcga_crc_clinical_v01.rds")

# ---- 1) Get SE (cache) ----

se <- NULL
if (file.exists(se_rds)) {
  log_message("Loading cached SE from: ", se_rds)
  se <- readRDS(se_rds)
} else {
  log_message("Querying TCGA COAD+READ (STAR - Counts)...")
  
  query <- GDCquery(
    project = c("TCGA-COAD", "TCGA-READ"),
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  GDCdownload(query)
  se <- tryCatch(
    { GDCprepare(query) },
    error = function(e) {
      log_message("ERROR in GDCprepare: ", e$message)
      stop("GDCprepare failed.")
    }
  )
  
  saveRDS(se, se_rds)
  log_message("SE saved to: ", se_rds)
}

if (!inherits(se, "SummarizedExperiment")) {
  stop("Object 'se' is not a SummarizedExperiment.")
}

assays_available <- assayNames(se)
expr_assay_name <- if ("unstranded" %in% assays_available) "unstranded" else assays_available[1]
log_message("Assays: ", paste(assays_available, collapse = ", "), " | Using: ", expr_assay_name)

# ---- 2) Clinical + survival ----

coldata <- as.data.frame(colData(se))

# Primary tumor
if ("shortLetterCode" %in% colnames(coldata)) {
  keep_primary <- coldata$shortLetterCode == "TP"
} else {
  keep_primary <- rep(TRUE, nrow(coldata))
}

se      <- se[, keep_primary]
coldata <- coldata[keep_primary, , drop = FALSE]
log_message("After primary tumor filter: ", ncol(se), " samples.")

# OS time: prefer days_to_death; if NA, use days_to_last_follow_up
dtd  <- suppressWarnings(as.numeric(coldata$days_to_death))
dtlf <- suppressWarnings(as.numeric(coldata$days_to_last_follow_up))
OS_time <- ifelse(!is.na(dtd), dtd, dtlf)

# Event
if ("vital_status" %in% colnames(coldata)) {
  OS_event <- ifelse(coldata$vital_status == "Dead", 1L, 0L)
} else {
  stop("vital_status not found for OS_event.")
}

keep_surv <- !is.na(OS_time) & !is.na(OS_event) & OS_time > 0
se        <- se[, keep_surv]
coldata   <- coldata[keep_surv, , drop = FALSE]
OS_time   <- OS_time[keep_surv]
OS_event  <- OS_event[keep_surv]

log_message("After OS filtering: ", length(OS_time), " samples.")

# Stage simplification
stage_candidates <- c(
  "ajcc_pathologic_tumor_stage",
  "pathologic_stage",
  "ajcc_pathologic_stage"
)
stage_col_present <- stage_candidates[stage_candidates %in% colnames(coldata)]
if (length(stage_col_present) > 0) {
  stage_raw <- toupper(as.character(coldata[[stage_col_present[1]]]))
  stage_simplified <- case_when(
    grepl("III", stage_raw) ~ "III/IV",
    grepl("IV", stage_raw) ~ "III/IV",
    grepl("I", stage_raw) & !grepl("III|IV", stage_raw) ~ "I/II",
    grepl("II", stage_raw) & !grepl("III|IV", stage_raw) ~ "I/II",
    TRUE ~ NA_character_
  )
} else {
  stage_simplified <- rep(NA_character_, nrow(coldata))
}

coldata$OS_time          <- OS_time
coldata$OS_event         <- OS_event
coldata$stage_simplified <- stage_simplified

log_message("Stage simplified available for ", sum(!is.na(stage_simplified)), " samples.")

# ---- 3) Expression → gene symbol level ----

expr_counts <- assay(se, expr_assay_name)
row_anno <- as.data.frame(rowData(se))

if ("external_gene_name" %in% colnames(row_anno)) {
  gene_symbols_all <- as.character(row_anno$external_gene_name)
} else if ("gene_name" %in% colnames(row_anno)) {
  gene_symbols_all <- as.character(row_anno$gene_name)
} else {
  log_message("WARNING: no gene symbol column; using rownames.")
  gene_symbols_all <- rownames(expr_counts)
}

dge <- DGEList(counts = expr_counts)

keep_genes <- filterByExpr(dge, group = OS_event)
dge <- dge[keep_genes, ]
gene_symbols <- gene_symbols_all[keep_genes]

dge <- calcNormFactors(dge)
design <- model.matrix(~ 1, data = coldata)
voom_res <- voom(dge, design = design, plot = FALSE)
expr_norm <- voom_res$E

valid <- !is.na(gene_symbols) & gene_symbols != ""
expr_norm <- expr_norm[valid, , drop = FALSE]
gene_symbols <- gene_symbols[valid]

tbl <- table(gene_symbols)
expr_sum <- rowsum(expr_norm, group = gene_symbols)
expr_norm_symbol <- sweep(expr_sum, 1, as.numeric(tbl[rownames(expr_sum)]), "/")

expr_norm <- expr_norm_symbol
log_message("Expression at gene-symbol level: ", nrow(expr_norm), " genes x ", ncol(expr_norm), " samples.")

# ---- 4) Save outputs ----

saveRDS(expr_norm, expr_rds)
saveRDS(coldata,   clin_rds)

log_message("Saved expr_norm to: ", expr_rds)
log_message("Saved clinical to: ", clin_rds)
log_message("03_TCGA_download_preprocess.R finished.")
