## 30_make_Figure3_clinical_axes_final_v03.R
## Final corrected multi-panel Figure 3
## Panel D uses recomputed Supplementary Table S1

rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

cat("[", Sys.time(), "] Figure 3 v03 script started.\n", sep = "")

## ---------------------------------------------------------
## 1) Paths
## ---------------------------------------------------------

project_root <- "."
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

dir_in  <- file.path(project_root, "results", "figure_build_inputs")
dir_tab <- file.path(project_root, "results", "tables")
dir_out <- file.path(project_root, "results", "figures_main")

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

path_fig3A <- file.path(dir_in,  "fig3A_stage_delta_input_v01.csv")
path_fig3B <- file.path(dir_in,  "fig3B_SCFA_stage_long_v01.csv")
path_fig3C <- file.path(dir_in,  "fig3C_SCFA_OS_forest_input_v01.csv")
path_s1    <- file.path(dir_tab, "Supplementary_Table_S1_multivariable_Cox_recomputed_v02.csv")

needed <- c(path_fig3A, path_fig3B, path_fig3C, path_s1)
missing <- needed[!file.exists(needed)]
if (length(missing) > 0) {
  stop("Missing required input file(s):\n", paste(missing, collapse = "\n"))
}

## ---------------------------------------------------------
## 2) Helpers
## ---------------------------------------------------------

axis_levels <- c("SCFA", "Tryptophan", "Polyamine", "Ferroptosis")
cohort_levels <- c("TCGA", "GSE39582")

cohort_cols <- c(
  "TCGA" = "#00A7B3",
  "GSE39582" = "#F26B61"
)

stage_cols <- c(
  "I/II" = "#F26B61",
  "III/IV" = "#00A7B3"
)

axis_clean <- function(x) {
  dplyr::recode(
    x,
    "SCFA_axis" = "SCFA",
    "Trp_axis" = "Tryptophan",
    "Polyamine_axis" = "Polyamine",
    "Ferroptosis_axis" = "Ferroptosis",
    "Tryptophan / AhR" = "Tryptophan",
    "Ferroptosis / redox" = "Ferroptosis",
    .default = x
  )
}

format_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

sig_label <- function(fdr) {
  dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 0.001 ~ "***",
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  )
}

parse_p <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  out <- suppressWarnings(as.numeric(x))
  lt <- str_detect(x, "^<")
  out[lt] <- suppressWarnings(as.numeric(str_remove(x[lt], "^<"))) / 2
  out
}

parse_hr_ci <- function(x) {
  m <- str_match(
    as.character(x),
    "^\\s*([0-9.]+)\\s*\\(([0-9.]+)\\s+to\\s+([0-9.]+)\\)\\s*$"
  )
  tibble(
    HR = as.numeric(m[, 2]),
    lower_CI = as.numeric(m[, 3]),
    upper_CI = as.numeric(m[, 4])
  )
}

theme_postbio <- function(base_size = 9.5) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, color = "grey15"),
      strip.background = element_rect(fill = "grey94", color = "grey35", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.30),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.20),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
}

## ---------------------------------------------------------
## 3) Read inputs
## ---------------------------------------------------------

fig3A <- read_csv(path_fig3A, show_col_types = FALSE)
fig3B <- read_csv(path_fig3B, show_col_types = FALSE)
fig3C <- read_csv(path_fig3C, show_col_types = FALSE)
s1    <- read_csv(path_s1, show_col_types = FALSE)

cat("[", Sys.time(), "] Inputs loaded.\n", sep = "")

## ---------------------------------------------------------
## 4) Clean inputs
## ---------------------------------------------------------

fig3A <- fig3A %>%
  mutate(
    axis_label = factor(axis_clean(axis), levels = axis_levels),
    cohort = factor(cohort, levels = cohort_levels),
    sig = sig_label(FDR_stage),
    is_sig = FDR_stage < 0.05,
    vjust_label = ifelse(delta_median < -0.17, 1.55, -0.85)
  )

fig3B <- fig3B %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    stage_simplified = factor(stage_simplified, levels = c("I/II", "III/IV"))
  ) %>%
  filter(!is.na(score), !is.na(stage_simplified), !is.na(cohort))

fig3C <- fig3C %>%
  mutate(
    cohort = factor(cohort, levels = rev(cohort_levels)),
    label = paste0(
      "HR ", sprintf("%.2f", HR),
      " (", sprintf("%.2f", lower_CI),
      " to ", sprintf("%.2f", upper_CI),
      "); FDR ", format_p(FDR_OS)
    )
  )

s1_parsed <- bind_cols(
  s1,
  parse_hr_ci(s1$`Adjusted HR (95% CI)`)
)

fig3D <- s1_parsed %>%
  mutate(
    cohort = factor(Cohort, levels = cohort_levels),
    axis_label = factor(axis_clean(Axis), levels = rev(axis_levels)),
    wald_p_num = parse_p(`Wald p value`),
    is_sig = wald_p_num < 0.05
  )

scfa_stage_labels <- fig3A %>%
  filter(axis_label == "SCFA") %>%
  mutate(
    label = paste0("FDR ", format_p(FDR_stage)),
    y = 1.62
  )

## ---------------------------------------------------------
## 5) Panel A
## ---------------------------------------------------------

pA <- ggplot(fig3A, aes(x = axis_label, y = delta_median, color = cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45, color = "black") +
  geom_point(
    aes(alpha = is_sig),
    size = 3.6,
    position = position_dodge(width = 0.45)
  ) +
  geom_text(
    aes(label = sig, group = cohort, vjust = vjust_label),
    position = position_dodge(width = 0.45),
    size = 3.0,
    color = "black",
    show.legend = FALSE
  ) +
  scale_color_manual(values = cohort_cols, name = "Cohort") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.30), guide = "none") +
  scale_y_continuous(
    limits = c(-0.235, 0.035),
    breaks = c(-0.20, -0.15, -0.10, -0.05, 0)
  ) +
  labs(
    title = "A",
    x = NULL,
    y = "Median difference\n(stage III/IV minus I/II)"
  ) +
  theme_postbio(base_size = 9.5) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 9.2)
  )

## ---------------------------------------------------------
## 6) Panel B
## ---------------------------------------------------------

pB <- ggplot(fig3B, aes(x = stage_simplified, y = score, fill = stage_simplified)) +
  geom_violin(
    width = 0.78,
    trim = TRUE,
    alpha = 0.82,
    linewidth = 0.32,
    color = "grey25"
  ) +
  geom_boxplot(
    width = 0.16,
    outlier.shape = NA,
    alpha = 0.85,
    linewidth = 0.32,
    color = "grey15"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.40, color = "black") +
  facet_wrap(~ cohort, nrow = 1) +
  geom_text(
    data = scfa_stage_labels,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.8,
    color = "grey10"
  ) +
  scale_fill_manual(values = stage_cols, name = "Stage") +
  labs(
    title = "B",
    x = "Stage group",
    y = "SCFA axis score"
  ) +
  coord_cartesian(ylim = c(-2.55, 1.75)) +
  theme_postbio(base_size = 9.5) +
  theme(
    legend.position = "none",
    panel.spacing = unit(0.45, "lines")
  )

## ---------------------------------------------------------
## 7) Panel C
## ---------------------------------------------------------

pC <- ggplot(fig3C, aes(x = HR, y = cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, color = "black") +
  geom_segment(
    aes(x = lower_CI, xend = upper_CI, y = cohort, yend = cohort),
    linewidth = 0.75,
    color = "black"
  ) +
  geom_point(size = 3.6, color = "black") +
  geom_text(
    aes(label = label),
    nudge_y = -0.16,
    hjust = 0.5,
    size = 2.45,
    color = "grey15"
  ) +
  scale_x_log10(
    limits = c(0.25, 1.35),
    breaks = c(0.25, 0.50, 0.75, 1.00, 1.25),
    labels = c("0.25", "0.50", "0.75", "1.00", "1.25")
  ) +
  labs(
    title = "C",
    x = "Hazard ratio for OS\n(per unit SCFA axis)",
    y = NULL
  ) +
  theme_postbio(base_size = 9.5) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

## ---------------------------------------------------------
## 8) Panel D
## ---------------------------------------------------------

pD <- ggplot(fig3D, aes(x = HR, y = axis_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, color = "black") +
  geom_segment(
    aes(x = lower_CI, xend = upper_CI, y = axis_label, yend = axis_label),
    linewidth = 0.70,
    color = "black"
  ) +
  geom_point(
    aes(fill = is_sig),
    shape = 21,
    size = 3.2,
    color = "black",
    stroke = 0.35
  ) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "white"),
    labels = c("FALSE" = "p >= 0.05", "TRUE" = "p < 0.05"),
    name = "Axis term"
  ) +
  scale_x_log10(
    limits = c(0.35, 2.8),
    breaks = c(0.5, 1.0, 2.0),
    labels = c("0.5", "1.0", "2.0")
  ) +
  labs(
    title = "D",
    x = "Adjusted hazard ratio for OS",
    y = NULL
  ) +
  theme_postbio(base_size = 9.5) +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(-4, 0, 0, 0),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.45, "lines")
  )

## ---------------------------------------------------------
## 9) Compose and save
## ---------------------------------------------------------

fig3 <- (pA / pB / (pC | pD)) +
  plot_layout(
    heights = c(0.78, 1.05, 1.05),
    widths = c(1, 1)
  )

out_pdf <- file.path(dir_out, "Figure3_clinical_axes_final_v03.pdf")
out_png <- file.path(dir_out, "Figure3_clinical_axes_final_v03.png")

ggsave(
  filename = out_pdf,
  plot = fig3,
  width = 11.8,
  height = 9.4,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = out_png,
  plot = fig3,
  width = 11.8,
  height = 9.4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(file.path(dir_out, "Figure3A_stage_delta_final_v03.pdf"), pA, width = 6.8, height = 3.1, device = cairo_pdf)
ggsave(file.path(dir_out, "Figure3B_SCFA_stage_final_v03.pdf"), pB, width = 7.2, height = 3.3, device = cairo_pdf)
ggsave(file.path(dir_out, "Figure3C_SCFA_OS_univ_final_v03.pdf"), pC, width = 4.8, height = 3.2, device = cairo_pdf)
ggsave(file.path(dir_out, "Figure3D_adjusted_OS_axes_final_v03.pdf"), pD, width = 6.3, height = 3.2, device = cairo_pdf)

cat("[", Sys.time(), "] Figure 3 v03 saved:\n", sep = "")
cat("PDF: ", out_pdf, "\n", sep = "")
cat("PNG: ", out_png, "\n", sep = "")
cat("[", Sys.time(), "] Figure 3 v03 script finished.\n", sep = "")
