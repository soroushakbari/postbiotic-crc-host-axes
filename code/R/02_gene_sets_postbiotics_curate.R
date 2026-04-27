## 02_gene_sets_postbiotics_curate.R
## Curate postbiotic axes and gene sets into a clean R list

source("./code/R/_config_Oncobiome.R")

# Packages
cran_pkgs <- c("readxl")
bioc_pkgs <- c("org.Hs.eg.db", "AnnotationDbi")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)

log_file <- file.path(dir_logs, paste0("02_gene_sets_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_message <- function(...) {
  msg <- paste0("[", Sys.time(), "] ", paste0(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_message("02_gene_sets_postbiotics_curate.R started.")

# ---- 1) Paths ----

axes_xlsx  <- file.path(dir_data_raw, "annotations", "postbiotic_axes_v01.xlsx")
genes_csv  <- file.path(dir_data_raw, "annotations", "postbiotic_gene_sets_v01.csv")
out_rds    <- file.path(dir_data_processed, "annotations", "postbiotic_gene_sets_list_v01.rds")
summary_csv <- file.path(dir_results, "tables", "postbiotic_gene_sets_summary_v01.csv")
dir.create(dirname(summary_csv), recursive = TRUE, showWarnings = FALSE)

# ---- 2) Read annotation files ----

if (!file.exists(axes_xlsx)) {
  stop("Axes file not found: ", axes_xlsx)
}
if (!file.exists(genes_csv)) {
  stop("Gene sets file not found: ", genes_csv)
}

axes_df <- readxl::read_excel(axes_xlsx, sheet = "axes")
axes_df <- as.data.frame(axes_df, stringsAsFactors = FALSE)

genes_df <- read.csv(genes_csv, stringsAsFactors = FALSE)

log_message("Read ", nrow(axes_df), " axes and ", nrow(genes_df), " gene annotations.")

# ---- 3) Basic sanity checks ----

required_axes_cols <- c("axis_id", "class", "short_description")
missing_axes <- setdiff(required_axes_cols, colnames(axes_df))
if (length(missing_axes) > 0) {
  stop("Missing required columns in axes sheet: ", paste(missing_axes, collapse = ", "))
}

required_gene_cols <- c("axis_id", "gene_symbol", "evidence_type", "priority")
missing_gene <- setdiff(required_gene_cols, colnames(genes_df))
if (length(missing_gene) > 0) {
  stop("Missing required columns in gene sets file: ", paste(missing_gene, collapse = ", "))
}

axes_df$axis_id <- trimws(as.character(axes_df$axis_id))
genes_df$axis_id <- trimws(as.character(genes_df$axis_id))
genes_df$gene_symbol <- toupper(trimws(as.character(genes_df$gene_symbol)))

# Axis IDs consistency
valid_axes <- unique(axes_df$axis_id)
bad_axes <- setdiff(unique(genes_df$axis_id), valid_axes)
if (length(bad_axes) > 0) {
  log_message("WARNING: gene_sets contain axis_id not in axes file: ",
              paste(bad_axes, collapse = ", "))
}

# ---- 4) Map gene symbols to Entrez / check validity ----

sym <- genes_df$gene_symbol
mapped <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(sym),
  columns = c("ENTREZID"),
  keytype = "SYMBOL"
)

mapped <- mapped[!is.na(mapped$ENTREZID), , drop = FALSE]

valid_symbols <- unique(mapped$SYMBOL)
genes_df$valid_symbol <- genes_df$gene_symbol %in% valid_symbols

n_invalid <- sum(!genes_df$valid_symbol)
if (n_invalid > 0) {
  log_message("WARNING: ", n_invalid, " gene symbols did not map to org.Hs.eg.db. They will be dropped.")
}

genes_df_valid <- genes_df[genes_df$valid_symbol, ]

# ---- 5) Build gene_sets list ----

gene_sets <- split(genes_df_valid$gene_symbol, genes_df_valid$axis_id)
gene_sets <- lapply(gene_sets, function(x) sort(unique(x)))

# Drop very small sets (less than 5 genes)
sizes <- sapply(gene_sets, length)
if (any(sizes < 5)) {
  log_message("Dropping axes with <5 genes: ",
              paste(names(sizes)[sizes < 5], collapse = ", "))
  gene_sets <- gene_sets[sizes >= 5]
  sizes <- sapply(gene_sets, length)
}

# ---- 6) Summary table ----

summary_df <- data.frame(
  axis_id   = names(gene_sets),
  n_genes   = as.integer(sizes),
  class     = axes_df$class[match(names(gene_sets), axes_df$axis_id)],
  short_description = axes_df$short_description[match(names(gene_sets), axes_df$axis_id)],
  stringsAsFactors = FALSE
)

write.csv(summary_df, summary_csv, row.names = FALSE)
saveRDS(gene_sets, out_rds)

log_message("Curated gene_sets saved to: ", out_rds)
log_message("Summary written to: ", summary_csv)
log_message("02_gene_sets_postbiotics_curate.R finished.")
