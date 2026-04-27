setwd(".")

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(readr)
  library(tibble)
})

fmt_num <- function(x, digits = 2) sprintf(paste0("%.", digits, "f"), x)

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

pretty_axis <- function(x) {
  dplyr::recode(
    x,
    "SCFA_axis" = "SCFA",
    "Trp_axis" = "Tryptophan / AhR",
    "Polyamine_axis" = "Polyamine",
    "Ferroptosis_axis" = "Ferroptosis / redox",
    .default = x
  )
}

fit_axis_cox <- function(clin, scores, cohort_name, covariates_label, covariates_formula) {
  
  common <- intersect(rownames(clin), rownames(scores))
  
  dat <- cbind(
    clin[common, , drop = FALSE],
    as.data.frame(scores[common, , drop = FALSE])
  )
  
  dat$stage_simplified <- factor(dat$stage_simplified, levels = c("I/II", "III/IV"))
  
  if ("age_at_diagnosis" %in% colnames(dat)) {
    dat$age_years <- as.numeric(dat$age_at_diagnosis) / 365.25
  }
  
  axes <- c("SCFA_axis", "Trp_axis", "Polyamine_axis", "Ferroptosis_axis")
  
  bind_rows(lapply(axes, function(ax) {
    
    dat_ax <- dat %>%
      filter(
        !is.na(OS_time),
        !is.na(OS_event),
        OS_time > 0,
        !is.na(.data[[ax]])
      )
    
    for (cc in covariates_formula) {
      dat_ax <- dat_ax %>% filter(!is.na(.data[[cc]]))
    }
    
    dat_ax$axis_score <- as.numeric(dat_ax[[ax]])
    
    rhs <- paste(c("axis_score", covariates_formula), collapse = " + ")
    form <- as.formula(paste0("Surv(OS_time, OS_event) ~ ", rhs))
    
    fit <- coxph(form, data = dat_ax)
    sm <- summary(fit)
    
    hr <- sm$coefficients["axis_score", "exp(coef)"]
    p  <- sm$coefficients["axis_score", "Pr(>|z|)"]
    lo <- sm$conf.int["axis_score", "lower .95"]
    hi <- sm$conf.int["axis_score", "upper .95"]
    
    tibble(
      cohort = cohort_name,
      axis = pretty_axis(ax),
      complete_case_n = nrow(dat_ax),
      adjustment_covariates = covariates_label,
      hr = hr,
      lower_ci = lo,
      upper_ci = hi,
      adjusted_hr_95_ci = paste0(
        fmt_num(hr, 2), " (", fmt_num(lo, 2), " to ", fmt_num(hi, 2), ")"
      ),
      wald_p_value = fmt_p(p)
    )
  }))
}

clin_tcga <- readRDS("data/processed/tcga_clinical/tcga_crc_clinical_v01.rds")
scores_tcga <- readRDS("data/processed/tcga_scores/tcga_crc_postbiotic_pathway_scores_v01.rds")

clin_gse <- readRDS("data/processed/external/GSE39582/GSE39582_clinical_v02.rds")
scores_gse <- readRDS("data/processed/external/GSE39582/GSE39582_postbiotic_pathway_scores_v02.rds")

s1_tcga <- fit_axis_cox(
  clin = clin_tcga,
  scores = scores_tcga,
  cohort_name = "TCGA",
  covariates_label = "Stage + age at diagnosis",
  covariates_formula = c("stage_simplified", "age_years")
)

s1_gse <- fit_axis_cox(
  clin = clin_gse,
  scores = scores_gse,
  cohort_name = "GSE39582",
  covariates_label = "Stage",
  covariates_formula = c("stage_simplified")
)

s1_final <- bind_rows(s1_tcga, s1_gse) %>%
  select(
    cohort,
    axis,
    complete_case_n,
    adjustment_covariates,
    adjusted_hr_95_ci,
    wald_p_value
  ) %>%
  rename(
    Cohort = cohort,
    Axis = axis,
    `Complete-case n` = complete_case_n,
    `Adjustment covariates` = adjustment_covariates,
    `Adjusted HR (95% CI)` = adjusted_hr_95_ci,
    `Wald p value` = wald_p_value
  )

print(s1_final, n = Inf)

out_csv <- "results/tables/Supplementary_Table_S1_multivariable_Cox_recomputed_v02.csv"
write_csv(s1_final, out_csv)

cat("\nWritten to:", out_csv, "\n")
