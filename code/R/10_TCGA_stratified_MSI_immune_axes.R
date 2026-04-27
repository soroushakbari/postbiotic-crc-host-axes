## 10_TCGA_stratified_MSI_immune_axes.R
## هدف: 
##  - مقایسه‌ی محورهای postbiotic بین MSI-H و MSS/MSI-L
##  - مقایسه‌ی محورهای postbiotic بین immune-hot و immune-cold (براساس Cytolytic_activity)
## خروجی: دو تا CSV در results/tables

rm(list = ls())

message("[", Sys.time(), "] 10_TCGA_stratified_MSI_immune_axes.R started.")

## 0) config و پکیج‌ها ----
project_root <- "."
config_path  <- file.path(project_root, "code", "R", "_config_Oncobiome.R")
if (!file.exists(config_path)) {
  stop("Config file not found at: ", config_path)
}
source(config_path)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

dir_clin    <- file.path(project_root, "data", "processed", "tcga_clinical")
dir_scores  <- file.path(project_root, "data", "processed", "tcga_scores")
dir_results <- file.path(project_root, "results", "tables")

if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE, showWarnings = FALSE)

## 1) لود کردن clinical, axes, TME ----

## clinical
clin_tcga_path <- file.path(dir_clin, "tcga_crc_clinical_v01.rds")
if (!file.exists(clin_tcga_path)) {
  stop("TCGA clinical file not found: ", clin_tcga_path)
}
clin_tcga <- readRDS(clin_tcga_path)

## axis scores
axes_tcga_path <- file.path(dir_scores, "tcga_crc_postbiotic_pathway_scores_v01.rds")
if (!file.exists(axes_tcga_path)) {
  stop("TCGA axis scores file not found: ", axes_tcga_path)
}
axes_tcga <- readRDS(axes_tcga_path)

## TME signatures
tme_path <- file.path(dir_scores, "tcga_crc_TME_signatures_v01.rds")
if (!file.exists(tme_path)) {
  stop("TCGA TME signatures file not found: ", tme_path)
}
tme_tcga <- readRDS(tme_path)

message("[", Sys.time(), "] Loaded TCGA clinical, axis scores, and TME signatures.")

## 2) align sample IDs ----

## فرض: rownames(clin_tcga) و rownames(axes_tcga) و rownames(tme_tcga) = sample barcodes
common_ids <- Reduce(
  intersect,
  list(rownames(clin_tcga), rownames(axes_tcga), rownames(tme_tcga))
)

if (length(common_ids) < 100) {
  stop("Too few overlapping samples across clinical, axes, and TME: n = ", length(common_ids))
}

clin   <- clin_tcga[common_ids, , drop = FALSE]
axes   <- axes_tcga[common_ids, , drop = FALSE]
tme    <- tme_tcga[common_ids, , drop = FALSE]

message("[", Sys.time(), "] Using ", length(common_ids), " overlapping TCGA samples.")

axes_to_use <- c("SCFA_axis", "Polyamine_axis", "Ferroptosis_axis", "Trp_axis")

## helper برای Wilcoxon و خلاصه
wilcox_summary <- function(values, group_factor, g1, g2) {
  v1 <- values[group_factor == g1]
  v2 <- values[group_factor == g2]
  if (length(v1) < 10 || length(v2) < 10) {
    return(list(
      n1 = length(v1),
      n2 = length(v2),
      med1 = NA_real_,
      med2 = NA_real_,
      diff_median = NA_real_,
      p_value = NA_real_
    ))
  }
  test <- suppressWarnings(wilcox.test(v1, v2, exact = FALSE))
  list(
    n1 = length(v1),
    n2 = length(v2),
    med1 = median(v1, na.rm = TRUE),
    med2 = median(v2, na.rm = TRUE),
    diff_median = median(v2, na.rm = TRUE) - median(v1, na.rm = TRUE),
    p_value = unname(test$p.value)
  )
}

## 3) MSI-H vs MSS/MSI-L ----

msi_res <- NULL

if ("paper_MSI_status" %in% colnames(clin)) {
  msi_raw <- as.character(clin$paper_MSI_status)
  
  msi_simpl <- ifelse(
    grepl("MSI-H", msi_raw, ignore.case = TRUE),
    "MSI-H",
    ifelse(is.na(msi_raw) | msi_raw == "", NA, "MSS_or_MSI-L")
  )
  
  table_msi <- table(msi_simpl, useNA = "ifany")
  message("[", Sys.time(), "] MSI simplified table: ")
  print(table_msi)
  
  msi_factor <- factor(msi_simpl, levels = c("MSS_or_MSI-L", "MSI-H"))
  
  ## فقط نمونه‌هایی که MSI تعریف‌شده دارند
  keep_msi <- !is.na(msi_factor)
  clin_msi <- clin[keep_msi, , drop = FALSE]
  axes_msi <- axes[keep_msi, , drop = FALSE]
  msi_factor <- msi_factor[keep_msi]
  
  message("[", Sys.time(), "] n with defined MSI status: ", length(msi_factor))
  
  res_list <- list()
  for (ax in axes_to_use) {
    if (!ax %in% colnames(axes_msi)) {
      warning("Axis ", ax, " not found in axes; skipping MSI analysis for this axis.")
      next
    }
    vals <- axes_msi[, ax]
    tmp <- wilcox_summary(vals, msi_factor, g1 = "MSS_or_MSI-L", g2 = "MSI-H")
    res_list[[ax]] <- data.frame(
      axis         = ax,
      group_ref    = "MSS_or_MSI-L",
      group_comp   = "MSI-H",
      n_ref        = tmp$n1,
      n_comp       = tmp$n2,
      median_ref   = tmp$med1,
      median_comp  = tmp$med2,
      diff_median  = tmp$diff_median,  ## MSI-H - MSS_or_MSI-L
      p_value      = tmp$p_value,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(res_list) > 0) {
    msi_res <- bind_rows(res_list) %>%
      mutate(
        FDR = p.adjust(p_value, method = "BH")
      )
    out_msi <- file.path(dir_results, "tcga_postbiotic_axes_MSI_assoc_v01.csv")
    write_csv(msi_res, out_msi)
    message("[", Sys.time(), "] MSI association results written to: ", out_msi)
  } else {
    warning("No MSI results produced.")
  }
  
} else {
  warning("paper_MSI_status column not found in clinical data; skipping MSI analysis.")
}

## 4) immune-hot vs immune-cold براساس Cytolytic_activity ----

immune_res <- NULL

if ("Cytolytic_activity" %in% colnames(tme)) {
  cyt <- tme[, "Cytolytic_activity"]
  if (!is.numeric(cyt)) {
    cyt <- as.numeric(cyt)
  }
  med_cyt <- median(cyt, na.rm = TRUE)
  immune_group <- ifelse(cyt > med_cyt, "Immune_high",
                         ifelse(is.na(cyt), NA, "Immune_low"))
  
  immune_factor <- factor(immune_group, levels = c("Immune_low", "Immune_high"))
  
  table_immune <- table(immune_factor, useNA = "ifany")
  message("[", Sys.time(), "] Immune group table (by cytolytic activity): ")
  print(table_immune)
  
  keep_imm <- !is.na(immune_factor)
  axes_imm <- axes[keep_imm, , drop = FALSE]
  immune_factor <- immune_factor[keep_imm]
  
  res_list2 <- list()
  for (ax in axes_to_use) {
    if (!ax %in% colnames(axes_imm)) {
      warning("Axis ", ax, " not found in axes; skipping immune-hot analysis for this axis.")
      next
    }
    vals <- axes_imm[, ax]
    tmp <- wilcox_summary(vals, immune_factor, g1 = "Immune_low", g2 = "Immune_high")
    res_list2[[ax]] <- data.frame(
      axis         = ax,
      group_ref    = "Immune_low",
      group_comp   = "Immune_high",
      n_ref        = tmp$n1,
      n_comp       = tmp$n2,
      median_ref   = tmp$med1,
      median_comp  = tmp$med2,
      diff_median  = tmp$diff_median,  ## Immune_high - Immune_low
      p_value      = tmp$p_value,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(res_list2) > 0) {
    immune_res <- bind_rows(res_list2) %>%
      mutate(
        FDR = p.adjust(p_value, method = "BH")
      )
    out_imm <- file.path(dir_results, "tcga_postbiotic_axes_immunehot_assoc_v01.csv")
    write_csv(immune_res, out_imm)
    message("[", Sys.time(), "] Immune-hot association results written to: ", out_imm)
  } else {
    warning("No immune-hot results produced.")
  }
  
} else {
  warning("Cytolytic_activity not found in TME signatures; skipping immune-hot analysis.")
}

message("[", Sys.time(), "] 10_TCGA_stratified_MSI_immune_axes.R finished.")
