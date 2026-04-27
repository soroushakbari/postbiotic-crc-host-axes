## 08b_GSE39582_scores_assoc_v02.R
## Compute postbiotic axis scores + stage/OS associations in GSE39582 (v02)

source("./code/R/_config_Oncobiome.R")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c("dplyr", "survival")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(dplyr)
library(survival)

log_file <- file.path(
  dir_logs,
  paste0("08b_GSE39582_scores_assoc_v02_",
         format(Sys.time(), "%Y%m%d_%H%M%S"),
         ".log")
)

log_message <- function(...) {
  msg <- paste0("[", Sys.time(), "] ", paste0(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_message("08b_GSE39582_scores_assoc_v02.R started.")

external_proc_dir <- file.path(dir_data_processed, "external", "GSE39582")
expr_rds_v02  <- file.path(external_proc_dir, "GSE39582_expr_symbol_log2_v02.rds")
clin_rds_v02  <- file.path(external_proc_dir, "GSE39582_clinical_v02.rds")
geneset_rds   <- file.path(dir_data_processed, "annotations",
                           "postbiotic_gene_sets_list_v01.rds")
scores_rds_v02 <- file.path(external_proc_dir,
                            "GSE39582_postbiotic_pathway_scores_v02.rds")

tables_dir <- file.path(dir_results, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
stage_tbl_v02 <- file.path(tables_dir,
                           "GSE39582_postbiotic_axes_stage_assoc_v02.csv")
os_tbl_v02    <- file.path(tables_dir,
                           "GSE39582_postbiotic_axes_OS_univ_v02.csv")

if (!file.exists(expr_rds_v02)) stop("Expression v02 RDS not found: ", expr_rds_v02)
if (!file.exists(clin_rds_v02)) stop("Clinical v02 RDS not found: ", clin_rds_v02)
if (!file.exists(geneset_rds)) stop("Gene sets RDS not found: ", geneset_rds)

expr_ext <- readRDS(expr_rds_v02)  # genes x samples (GSM IDs as colnames)
clin_ext <- readRDS(clin_rds_v02)  # samples as rownames
gene_sets <- readRDS(geneset_rds)

log_message(paste0("expr_ext dim: ",
                   nrow(expr_ext), " genes x ",
                   ncol(expr_ext), " samples."))
log_message(paste0("clin_ext dim: ",
                   nrow(clin_ext), " samples x ",
                   ncol(clin_ext), " columns."))
log_message("Example expr colnames: ",
            paste(head(colnames(expr_ext), 5), collapse = ", "))
log_message("Example clin rownames: ",
            paste(head(rownames(clin_ext), 5), collapse = ", "))

# -----------------------------------------------------------
# Z-score and compute scores
# -----------------------------------------------------------

log_message("Z-scoring expression (per gene)...")
z_expr <- t(scale(t(expr_ext)))

axes <- names(gene_sets)

score_mat <- sapply(axes, function(ax) {
  gs <- gene_sets[[ax]]
  gs <- intersect(gs, rownames(z_expr))
  if (length(gs) == 0) {
    log_message("WARNING: axis ", ax, " has 0 overlapping genes; scores will be NA.")
    return(rep(NA_real_, ncol(z_expr)))
  }
  colMeans(z_expr[gs, , drop = FALSE], na.rm = TRUE)
})

score_mat <- as.matrix(score_mat)
rownames(score_mat) <- colnames(z_expr)  # GSM IDs
colnames(score_mat) <- axes

log_message(paste0("score_mat dim: ",
                   nrow(score_mat), " samples x ",
                   ncol(score_mat), " axes."))
log_message("Example score rownames: ",
            paste(head(rownames(score_mat), 5), collapse = ", "))

saveRDS(score_mat, scores_rds_v02)
log_message("Saved scores to: ", scores_rds_v02)

# -----------------------------------------------------------
# Align clinical + scores
# -----------------------------------------------------------

common_samples <- intersect(rownames(clin_ext), rownames(score_mat))
log_message(paste0("Overlapping samples between clin_ext and scores: ",
                   length(common_samples)))

clin_sub   <- clin_ext[common_samples, , drop = FALSE]
scores_sub <- score_mat[common_samples, , drop = FALSE]

# -----------------------------------------------------------
# Stage association
# -----------------------------------------------------------

stage_mask <- !is.na(clin_sub$stage_simplified) &
  clin_sub$stage_simplified %in% c("I/II", "III/IV")

dat_stage    <- clin_sub[stage_mask, , drop = FALSE]
scores_stage <- scores_sub[stage_mask, , drop = FALSE]

log_message(paste0("Stage association samples: ", nrow(dat_stage)))

stage_res <- lapply(axes, function(ax) {
  x1 <- scores_stage[, ax][dat_stage$stage_simplified == "I/II"]
  x2 <- scores_stage[, ax][dat_stage$stage_simplified == "III/IV"]
  
  if (length(x1) > 5 && length(x2) > 5) {
    wt <- wilcox.test(x1, x2)
    delta <- median(x2, na.rm = TRUE) - median(x1, na.rm = TRUE)
    data.frame(
      axis = ax,
      n_I_II = sum(dat_stage$stage_simplified == "I/II"),
      n_III_IV = sum(dat_stage$stage_simplified == "III/IV"),
      median_I_II = median(x1, na.rm = TRUE),
      median_III_IV = median(x2, na.rm = TRUE),
      delta_median = delta,
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      axis = ax,
      n_I_II = sum(dat_stage$stage_simplified == "I/II"),
      n_III_IV = sum(dat_stage$stage_simplified == "III/IV"),
      median_I_II = NA_real_,
      median_III_IV = NA_real_,
      delta_median = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    )
  }
})

stage_res <- do.call(rbind, stage_res)
stage_res$FDR <- p.adjust(stage_res$p_value, method = "BH")
write.csv(stage_res, stage_tbl_v02, row.names = FALSE)
log_message("Wrote stage association table to: ", stage_tbl_v02)

# -----------------------------------------------------------
# OS – univariable Cox per axis
# -----------------------------------------------------------

OS_time  <- as.numeric(clin_sub$OS_time)
OS_event <- as.numeric(clin_sub$OS_event)

os_res <- lapply(axes, function(ax) {
  score_vec <- scores_sub[, ax]
  valid <- !is.na(score_vec) & !is.na(OS_time) & !is.na(OS_event) & OS_time > 0
  if (sum(valid) < 50) {
    return(data.frame(
      axis = ax,
      n = sum(valid),
      HR = NA_real_,
      lower_CI = NA_real_,
      upper_CI = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  fit <- coxph(Surv(OS_time[valid], OS_event[valid]) ~ score_vec[valid])
  s <- summary(fit)
  HR <- s$coef[,"exp(coef)"][1]
  CI <- s$conf.int[1, c("lower .95", "upper .95")]
  pval <- s$coef[,"Pr(>|z|)"][1]
  
  data.frame(
    axis = ax,
    n = sum(valid),
    HR = HR,
    lower_CI = CI[1],
    upper_CI = CI[2],
    p_value = pval,
    stringsAsFactors = FALSE
  )
})

os_res <- do.call(rbind, os_res)
os_res$FDR <- p.adjust(os_res$p_value, method = "BH")
write.csv(os_res, os_tbl_v02, row.names = FALSE)
log_message("Wrote OS association table to: ", os_tbl_v02)

log_message("08b_GSE39582_scores_assoc_v02.R finished.")
