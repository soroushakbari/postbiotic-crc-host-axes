## 04_TCGA_pathway_scores.R
## Compute postbiotic pathway scores (z-score-based) for TCGA CRC

source("./code/R/_config_Oncobiome.R")

cran_pkgs <- c("dplyr")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(dplyr)

log_file <- file.path(dir_logs, paste0("04_pathway_scores_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_message <- function(...) {
  msg <- paste0("[", Sys.time(), "] ", paste0(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_message("04_TCGA_pathway_scores.R started.")

expr_rds   <- file.path(dir_data_processed, "tcga_expression", "tcga_crc_expr_v01.rds")
geneset_rds <- file.path(dir_data_processed, "annotations", "postbiotic_gene_sets_list_v01.rds")
scores_rds <- file.path(dir_data_processed, "tcga_scores", "tcga_crc_postbiotic_pathway_scores_v01.rds")
summary_csv <- file.path(dir_results, "tables", "tcga_postbiotic_pathway_scores_summary_v01.csv")
dir.create(dirname(summary_csv), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expr_rds)) stop("Expression file not found: ", expr_rds)
if (!file.exists(geneset_rds)) stop("Gene sets file not found: ", geneset_rds)

expr_norm <- readRDS(expr_rds)   # genes x samples
gene_sets <- readRDS(geneset_rds)

log_message("Loaded expr_norm (", nrow(expr_norm), " genes x ", ncol(expr_norm),
            " samples) and ", length(gene_sets), " gene sets.")

# ---- 1) Z-score across samples ----

z_expr <- t(scale(t(expr_norm)))
log_message("Z-scoring done.")

# ---- 2) Compute scores (mean z per axis) ----

axes <- names(gene_sets)

score_mat <- sapply(axes, function(ax) {
  gs <- gene_sets[[ax]]
  gs <- intersect(gs, rownames(z_expr))
  if (length(gs) == 0) {
    log_message("WARNING: axis ", ax, " has 0 genes present in expr_norm.")
    return(rep(NA_real_, ncol(z_expr)))
  }
  colMeans(z_expr[gs, , drop = FALSE], na.rm = TRUE)
})

score_mat <- as.matrix(score_mat)
rownames(score_mat) <- colnames(expr_norm)
colnames(score_mat) <- axes

log_message("Score matrix dim: ", nrow(score_mat), " samples x ", ncol(score_mat), " axes.")

# ---- 3) Save outputs ----

saveRDS(score_mat, scores_rds)

summary_df <- data.frame(
  axis_id = axes,
  n_genes_in_set = sapply(gene_sets, length),
  n_genes_used   = sapply(axes, function(ax) length(intersect(gene_sets[[ax]], rownames(z_expr)))),
  stringsAsFactors = FALSE
)

write.csv(summary_df, summary_csv, row.names = FALSE)

log_message("Scores saved to: ", scores_rds)
log_message("Summary written to: ", summary_csv)
log_message("04_TCGA_pathway_scores.R finished.")
