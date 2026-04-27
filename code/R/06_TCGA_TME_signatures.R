## 06_TCGA_TME_signatures.R
## Build simple TME/immune signatures and correlate with postbiotic axes in TCGA CRC

source("./code/R/_config_Oncobiome.R")

# -----------------------------------------------------------
# Packages
# -----------------------------------------------------------

cran_pkgs <- c("dplyr", "ggplot2", "reshape2")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(dplyr)
library(ggplot2)
library(reshape2)

log_file <- file.path(dir_logs, paste0("06_TME_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_message <- function(...) {
  msg <- paste0("[", Sys.time(), "] ", paste0(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_message("06_TCGA_TME_signatures.R started.")

# -----------------------------------------------------------
# Paths
# -----------------------------------------------------------

expr_rds     <- file.path(dir_data_processed, "tcga_expression", "tcga_crc_expr_v01.rds")
clin_rds     <- file.path(dir_data_processed, "tcga_clinical",   "tcga_crc_clinical_v01.rds")
axes_scores_rds <- file.path(dir_data_processed, "tcga_scores",  "tcga_crc_postbiotic_pathway_scores_v01.rds")

tme_scores_rds  <- file.path(dir_data_processed, "tcga_scores",  "tcga_crc_TME_signatures_v01.rds")
tme_summary_csv <- file.path(dir_results, "tables", "tcga_TME_signatures_summary_v01.csv")
cor_csv         <- file.path(dir_results, "tables", "tcga_postbiotic_axes_TME_cor_v01.csv")

fig_dir <- file.path(dir_results, "figures_main")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(tme_summary_csv), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------
# Load data
# -----------------------------------------------------------

if (!file.exists(expr_rds)) stop("Expression file not found: ", expr_rds)
if (!file.exists(clin_rds)) stop("Clinical file not found: ", clin_rds)
if (!file.exists(axes_scores_rds)) stop("Axes scores file not found: ", axes_scores_rds)

expr_norm <- readRDS(expr_rds)       # genes x samples
clin      <- readRDS(clin_rds)       # data.frame, rownames = sample IDs
axes_scores <- readRDS(axes_scores_rds)  # samples x axes

log_message("Loaded expr_norm (", nrow(expr_norm), " genes x ", ncol(expr_norm),
            " samples); clinical (", nrow(clin), " samples); axes_scores (",
            nrow(axes_scores), " samples x ", ncol(axes_scores), " axes).")

# -----------------------------------------------------------
# Align samples
# -----------------------------------------------------------

common_samples <- Reduce(intersect, list(colnames(expr_norm), rownames(clin), rownames(axes_scores)))
if (length(common_samples) < 50) {
  stop("Too few overlapping samples: ", length(common_samples))
}

expr_norm   <- expr_norm[, common_samples, drop = FALSE]
clin        <- clin[common_samples, , drop = FALSE]
axes_scores <- axes_scores[common_samples, , drop = FALSE]

log_message("Using ", length(common_samples), " overlapping samples.")

# -----------------------------------------------------------
# Define simple TME / immune signatures (manual, canonical markers)
# -----------------------------------------------------------

tme_sets_manual <- list(
  CD8_T_cells = c("CD8A", "CD8B", "GZMB", "PRF1", "IFNG"),
  Tregs       = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"),
  Th1_like    = c("TBX21", "STAT1", "IFNG", "CXCR3"),
  Th17_like   = c("RORC", "IL17A", "IL17F", "CCR6"),
  NK_cells    = c("NKG7", "GNLY", "KLRD1", "KLRK1", "PRF1"),
  M1_macrophages = c("IRF5", "STAT1", "CD80", "CD86", "IL12B"),
  M2_macrophages = c("CD163", "MSR1", "MRC1", "IL10", "TGFB1"),
  Cytolytic_activity = c("GZMA", "GZMB", "PRF1", "GNLY", "NKG7"),
  Stromal_like = c("COL1A1", "COL3A1", "ACTA2", "TAGLN", "PDGFRB")
)

all_genes <- rownames(expr_norm)

tme_sets <- lapply(tme_sets_manual, function(gs) intersect(gs, all_genes))
tme_sets <- tme_sets[sapply(tme_sets, length) >= 3]

if (length(tme_sets) == 0) {
  stop("All TME gene sets are empty after intersecting with expr_norm.")
}

log_message("TME gene sets (non-empty):")
for (nm in names(tme_sets)) {
  log_message(" - ", nm, ": ", length(tme_sets[[nm]]), " genes")
}

tme_names <- names(tme_sets)

# -----------------------------------------------------------
# Compute TME scores: per-gene z-score across samples, mean z per set
# -----------------------------------------------------------

log_message("Computing z-score-based TME signatures...")

z_expr <- t(scale(t(expr_norm)))  # genes x samples

tme_score_mat <- sapply(tme_names, function(sig) {
  gs <- tme_sets[[sig]]
  gs <- intersect(gs, rownames(z_expr))
  if (length(gs) == 0) {
    return(rep(NA_real_, ncol(z_expr)))
  }
  colMeans(z_expr[gs, , drop = FALSE], na.rm = TRUE)
})

tme_score_mat <- as.matrix(tme_score_mat)
rownames(tme_score_mat) <- colnames(expr_norm)  # samples
colnames(tme_score_mat) <- tme_names

log_message("TME scores matrix dim: ", nrow(tme_score_mat), " samples x ",
            ncol(tme_score_mat), " signatures.")

# -----------------------------------------------------------
# Save TME scores and summary
# -----------------------------------------------------------

saveRDS(tme_score_mat, tme_scores_rds)

tme_summary <- data.frame(
  signature = tme_names,
  n_genes_defined = sapply(tme_sets_manual[tme_names], length),
  n_genes_used    = sapply(tme_sets, length),
  stringsAsFactors = FALSE
)

write.csv(tme_summary, tme_summary_csv, row.names = FALSE)
log_message("TME signatures summary written to: ", tme_summary_csv)
log_message("TME scores saved to: ", tme_scores_rds)

# -----------------------------------------------------------
# Correlation: postbiotic axes × TME signatures (Spearman)
# -----------------------------------------------------------

axes <- colnames(axes_scores)

# Build combined dataset
axes_df <- as.data.frame(axes_scores)
axes_df$sample <- rownames(axes_df)

tme_df <- as.data.frame(tme_score_mat)
tme_df$sample <- rownames(tme_df)

clin$sample <- rownames(clin)

dat <- clin %>%
  inner_join(axes_df, by = "sample") %>%
  inner_join(tme_df,  by = "sample")

log_message("Combined dataset for correlation: ", nrow(dat), " samples.")

cor_list <- list()

for (ax in axes) {
  for (tm in tme_names) {
    x <- dat[[ax]]
    y <- dat[[tm]]
    
    valid <- !is.na(x) & !is.na(y)
    if (sum(valid) < 30) {
      cor_list[[length(cor_list) + 1]] <- data.frame(
        axis = ax,
        TME_signature = tm,
        n = sum(valid),
        rho = NA_real_,
        p_value = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      ct <- suppressWarnings(cor.test(x[valid], y[valid], method = "spearman"))
      cor_list[[length(cor_list) + 1]] <- data.frame(
        axis = ax,
        TME_signature = tm,
        n = sum(valid),
        rho = unname(ct$estimate),
        p_value = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
}

cor_df <- do.call(rbind, cor_list)
cor_df$FDR <- p.adjust(cor_df$p_value, method = "BH")

write.csv(cor_df, cor_csv, row.names = FALSE)
log_message("Axes–TME correlation table written to: ", cor_csv)

# -----------------------------------------------------------
# Heatmap of rho (axes × TME)
# -----------------------------------------------------------

cor_mat <- matrix(NA_real_, nrow = length(axes), ncol = length(tme_names))
rownames(cor_mat) <- axes
colnames(cor_mat) <- tme_names

for (i in seq_len(nrow(cor_df))) {
  a  <- cor_df$axis[i]
  tm <- cor_df$TME_signature[i]
  r  <- cor_df$rho[i]
  if (!is.na(r)) {
    cor_mat[a, tm] <- r
  }
}

cor_plot_df <- reshape2::melt(cor_mat, varnames = c("axis", "TME_signature"), value.name = "rho")

p_heat <- ggplot(cor_plot_df, aes(x = TME_signature, y = axis, fill = rho)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey90"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Spearman correlation: postbiotic axes vs TME signatures (TCGA CRC)",
    x = "TME signature",
    y = "Postbiotic axis"
  )

heatmap_file <- file.path(fig_dir, "TCGA_postbiotic_axes_TME_correlation_heatmap_v01.pdf")
ggsave(heatmap_file, p_heat, width = 7, height = 4)
log_message("Correlation heatmap saved to: ", heatmap_file)

log_message("06_TCGA_TME_signatures.R finished.")
