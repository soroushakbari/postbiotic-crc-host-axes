## 05_TCGA_clinical_survival.R
## Clinical & OS associations for postbiotic pathway scores in TCGA CRC

source("./code/R/_config_Oncobiome.R")

cran_pkgs <- c("dplyr", "survival", "ggplot2")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(dplyr)
library(survival)
library(ggplot2)

log_file <- file.path(dir_logs, paste0("05_clinical_survival_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_message <- function(...) {
  msg <- paste0("[", Sys.time(), "] ", paste0(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_message("05_TCGA_clinical_survival.R started.")

# -----------------------------------------------------------
# Paths
# -----------------------------------------------------------

clin_rds   <- file.path(dir_data_processed, "tcga_clinical",  "tcga_crc_clinical_v01.rds")
scores_rds <- file.path(dir_data_processed, "tcga_scores",    "tcga_crc_postbiotic_pathway_scores_v01.rds")

if (!file.exists(clin_rds))   stop("Clinical file not found: ", clin_rds)
if (!file.exists(scores_rds)) stop("Scores file not found: ", scores_rds)

clin   <- readRDS(clin_rds)    # data.frame, rownames = sample IDs
scores <- readRDS(scores_rds)  # matrix: samples x axes

log_message("Loaded clinical (", nrow(clin), " samples) and scores (",
            nrow(scores), " samples x ", ncol(scores), " axes).")

# -----------------------------------------------------------
# Align samples
# -----------------------------------------------------------

common_samples <- intersect(rownames(clin), rownames(scores))
if (length(common_samples) < 50) {
  stop("Too few overlapping samples between clinical and scores: ", length(common_samples))
}

clin   <- clin[common_samples, , drop = FALSE]
scores <- scores[common_samples, , drop = FALSE]

axes <- colnames(scores)
log_message("Using ", length(common_samples), " samples and axes: ",
            paste(axes, collapse = ", "))

# Build combined dataset
scores_df <- as.data.frame(scores)
scores_df$sample <- rownames(scores_df)

clin$sample <- rownames(clin)
dat <- dplyr::left_join(clin, scores_df, by = "sample")

# -----------------------------------------------------------
# 1) Stage association (I/II vs III/IV)
# -----------------------------------------------------------

results_tables_dir <- file.path(dir_results, "tables")
dir.create(results_tables_dir, recursive = TRUE, showWarnings = FALSE)

stage_mask <- !is.na(dat$stage_simplified) & dat$stage_simplified %in% c("I/II", "III/IV")
dat_stage  <- dat[stage_mask, ]

log_message("Stage association: ", nrow(dat_stage),
            " samples with stage_simplified (I/II vs III/IV).")

stage_res <- lapply(axes, function(ax) {
  x1 <- dat_stage[[ax]][dat_stage$stage_simplified == "I/II"]
  x2 <- dat_stage[[ax]][dat_stage$stage_simplified == "III/IV"]
  
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
      median_I_II = NA,
      median_III_IV = NA,
      delta_median = NA,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    )
  }
})

stage_res <- do.call(rbind, stage_res)
stage_res$FDR <- p.adjust(stage_res$p_value, method = "BH")

stage_outfile <- file.path(results_tables_dir, "tcga_postbiotic_axes_stage_assoc_v01.csv")
write.csv(stage_res, stage_outfile, row.names = FALSE)
log_message("Stage association results written to: ", stage_outfile)

# -----------------------------------------------------------
# 2) OS – univariable Cox
# -----------------------------------------------------------

if (!("OS_time" %in% colnames(dat)) || !("OS_event" %in% colnames(dat))) {
  stop("OS_time / OS_event not found in clinical data.")
}

OS_time  <- as.numeric(dat$OS_time)
OS_event <- as.numeric(dat$OS_event)

surv_uni <- lapply(axes, function(ax) {
  score_vec <- dat[[ax]]
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

surv_uni <- do.call(rbind, surv_uni)
surv_uni$FDR <- p.adjust(surv_uni$p_value, method = "BH")

surv_uni_outfile <- file.path(results_tables_dir, "tcga_postbiotic_axes_OS_univ_v01.csv")
write.csv(surv_uni, surv_uni_outfile, row.names = FALSE)
log_message("Univariable OS Cox results written to: ", surv_uni_outfile)

# -----------------------------------------------------------
# 3) OS – multivariable Cox (score + age + stage)
# -----------------------------------------------------------

# Try to find an age column
age_candidates <- c("age_at_initial_pathologic_diagnosis", "age_at_diagnosis", "paper_age_at_diagnosis")
age_col_present <- age_candidates[age_candidates %in% colnames(dat)]

if (length(age_col_present) == 0) {
  log_message("No age column found; skipping multivariable Cox.")
} else {
  age_col <- age_col_present[1]
  log_message("Using age column for multivariable Cox: ", age_col)
  
  # Stage factor: I/II vs III/IV vs other/NA
  stage_mv <- factor(dat$stage_simplified, levels = c("I/II", "III/IV"))
  age_vec  <- suppressWarnings(as.numeric(dat[[age_col]]))
  
  surv_multi <- lapply(axes, function(ax) {
    score_vec <- dat[[ax]]
    valid <- !is.na(score_vec) & !is.na(OS_time) & !is.na(OS_event) &
      !is.na(age_vec) & !is.na(stage_mv) & OS_time > 0
    
    if (sum(valid) < 50) {
      return(data.frame(
        axis = ax,
        n = sum(valid),
        HR_score = NA_real_,
        lower_CI_score = NA_real_,
        upper_CI_score = NA_real_,
        p_value_score = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    fit <- coxph(
      Surv(OS_time[valid], OS_event[valid]) ~ score_vec[valid] + age_vec[valid] + stage_mv[valid]
    )
    s <- summary(fit)
    # coef for score is first term
    HR <- s$coef[1,"exp(coef)"]
    CI_low <- s$conf.int[1, "lower .95"]
    CI_up  <- s$conf.int[1, "upper .95"]
    pval <- s$coef[1,"Pr(>|z|)"]
    
    data.frame(
      axis = ax,
      n = sum(valid),
      HR_score = HR,
      lower_CI_score = CI_low,
      upper_CI_score = CI_up,
      p_value_score = pval,
      stringsAsFactors = FALSE
    )
  })
  
  surv_multi <- do.call(rbind, surv_multi)
  surv_multi$FDR_score <- p.adjust(surv_multi$p_value_score, method = "BH")
  
  surv_multi_outfile <- file.path(results_tables_dir, "tcga_postbiotic_axes_OS_multi_v01.csv")
  write.csv(surv_multi, surv_multi_outfile, row.names = FALSE)
  log_message("Multivariable OS Cox results written to: ", surv_multi_outfile)
}

# -----------------------------------------------------------
# 4) Simple stage boxplots (QC / for main figs later)
# -----------------------------------------------------------

fig_dir <- file.path(dir_results, "figures_main")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (nrow(dat_stage) > 0) {
  for (ax in axes) {
    if (!all(is.na(dat_stage[[ax]]))) {
      p <- ggplot(dat_stage, aes(x = stage_simplified, y = .data[[ax]])) +
        geom_boxplot() +
        labs(
          title = paste0("TCGA CRC: ", ax, " vs stage (I/II vs III/IV)"),
          x = "Stage group",
          y = "pathway score (mean gene z-score)"
        ) +
        theme_bw()
      
      ggsave(
        filename = file.path(fig_dir, paste0("TCGA_", ax, "_stage_boxplot_v01.pdf")),
        plot = p,
        width = 4, height = 4
      )
    }
  }
  log_message("Stage boxplots saved to: ", fig_dir)
}

log_message("05_TCGA_clinical_survival.R finished.")
