## 11_Fig2_literature_landscape_v02.R
## Polished Figure 2: literature evidence landscape of postbiotic-related axes in CRC

rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

cat("[", Sys.time(), "] Figure 2 v02 script started.\n", sep = "")

project_root <- "."
project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

infile <- file.path(
  project_root,
  "results", "tables",
  "Table3_postbiotic_evidence_matrix_for_Fig2_v01.csv"
)

out_dir <- file.path(project_root, "results", "figures_main")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(infile)) {
  stop("Input file not found: ", infile)
}

dat <- read_csv(infile, show_col_types = FALSE)

required_cols <- c(
  "study_id",
  "first_author_year",
  "axis",
  "study_type_simplified",
  "specimen_model_simplified",
  "endpoint_simplified",
  "evidence_strength",
  "include_in_fig2",
  "include_priority",
  "verification_flag",
  "short_label_for_plot"
)

missing_cols <- setdiff(required_cols, colnames(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

dat_fig2 <- dat %>%
  mutate(
    include_in_fig2 = as.character(include_in_fig2),
    evidence_strength = as.numeric(evidence_strength)
  ) %>%
  filter(tolower(include_in_fig2) %in% c("yes", "y", "true", "1"))

if (nrow(dat_fig2) == 0) {
  stop("No rows retained for Figure 2.")
}

axis_levels <- c(
  "SCFA",
  "Tryptophan-AhR",
  "Polyamine",
  "Ferroptosis-redox",
  "Mixed / multi-omics"
)

endpoint_levels <- c(
  "Diagnosis / early detection",
  "Risk / recurrence",
  "Stage / progression",
  "Prognosis / survival",
  "Therapy response",
  "TME / immune context",
  "Mechanistic evidence"
)

specimen_levels <- c(
  "Stool / faeces",
  "Blood / serum / urine",
  "Tumour / mucosa",
  "Preclinical model",
  "Multi-omics / literature"
)

axis_labels <- c(
  "SCFA" = "SCFA",
  "Tryptophan-AhR" = "Tryptophan / AhR",
  "Polyamine" = "Polyamine",
  "Ferroptosis-redox" = "Ferroptosis / redox",
  "Mixed / multi-omics" = "Mixed / multi-omics"
)

axis_labels_short <- c(
  "SCFA" = "SCFA",
  "Tryptophan-AhR" = "Tryptophan\nAhR",
  "Polyamine" = "Polyamine",
  "Ferroptosis-redox" = "Ferroptosis\nredox",
  "Mixed / multi-omics" = "Mixed\nmulti-omics"
)

endpoint_labels <- c(
  "Diagnosis / early detection" = "Diagnosis /\nearly detection",
  "Risk / recurrence" = "Risk /\nrecurrence",
  "Stage / progression" = "Stage /\nprogression",
  "Prognosis / survival" = "Prognosis /\nsurvival",
  "Therapy response" = "Therapy\nresponse",
  "TME / immune context" = "TME /\nimmune context",
  "Mechanistic evidence" = "Mechanistic\nevidence"
)

specimen_labels <- c(
  "Stool / faeces" = "Stool / faeces",
  "Blood / serum / urine" = "Blood / serum / urine",
  "Tumour / mucosa" = "Tumour / mucosa",
  "Preclinical model" = "Preclinical model",
  "Multi-omics / literature" = "Multi-omics / literature"
)

strength_labels <- c(
  "1" = "Review-level",
  "2" = "Preclinical",
  "3" = "Human observational",
  "4" = "Translational or multi-omics",
  "5" = "Clinical or interventional"
)

dat_fig2 <- dat_fig2 %>%
  mutate(
    axis = factor(axis, levels = axis_levels),
    endpoint_simplified = factor(endpoint_simplified, levels = endpoint_levels),
    specimen_model_simplified = factor(specimen_model_simplified, levels = specimen_levels),
    evidence_strength_factor = factor(evidence_strength, levels = 1:5)
  ) %>%
  filter(!is.na(axis))

theme_postbio_fig2 <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, color = "grey15"),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.30),
      panel.grid.minor = element_blank(),
      legend.title = element_text(face = "bold", size = base_size - 0.5),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.40, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

## Panel A
panelA_grid <- expand.grid(
  axis = factor(axis_levels, levels = axis_levels),
  endpoint_simplified = factor(endpoint_levels, levels = endpoint_levels)
) %>%
  as_tibble()

panelA_nonzero <- dat_fig2 %>%
  group_by(axis, endpoint_simplified) %>%
  summarise(
    n_studies = n(),
    mean_strength = mean(evidence_strength, na.rm = TRUE),
    .groups = "drop"
  )

panelA_df <- panelA_grid %>%
  left_join(panelA_nonzero, by = c("axis", "endpoint_simplified")) %>%
  mutate(
    n_studies = ifelse(is.na(n_studies), 0, n_studies),
    axis = factor(axis, levels = rev(axis_levels)),
    endpoint_simplified = factor(endpoint_simplified, levels = endpoint_levels)
  )

pA <- ggplot(panelA_df, aes(x = endpoint_simplified, y = axis)) +
  geom_tile(fill = "grey98", color = "grey88", linewidth = 0.35) +
  geom_point(
    data = filter(panelA_df, n_studies > 0),
    aes(size = n_studies, fill = mean_strength),
    shape = 21,
    color = "black",
    stroke = 0.35,
    alpha = 0.95
  ) +
  geom_text(
    data = filter(panelA_df, n_studies > 0),
    aes(label = n_studies),
    size = 2.9,
    color = "black"
  ) +
  scale_x_discrete(labels = endpoint_labels, drop = FALSE) +
  scale_y_discrete(labels = axis_labels, drop = FALSE) +
  scale_size_continuous(
    range = c(4, 12),
    breaks = c(1, 2, 3),
    name = "No. of\nstudies"
  ) +
  scale_fill_gradient(
    low = "#deebf7",
    high = "#08519c",
    limits = c(1, 5),
    breaks = 1:5,
    name = "Mean evidence\nstrength"
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  ) +
  labs(
    title = "Evidence density across axes and clinical-translational endpoints",
    x = NULL,
    y = NULL
  ) +
  theme_postbio_fig2(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "right",
    panel.grid = element_blank()
  )

## Panel B
panelB_df <- dat_fig2 %>%
  count(axis, specimen_model_simplified, name = "n") %>%
  group_by(axis) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  mutate(axis = factor(axis, levels = rev(axis_levels)))

pB <- ggplot(panelB_df, aes(x = frac, y = axis, fill = specimen_model_simplified)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.72) +
  geom_text(
    aes(label = ifelse(frac >= 0.12, n, "")),
    position = position_stack(vjust = 0.5),
    size = 2.8,
    color = "black"
  ) +
  scale_y_discrete(labels = axis_labels_short, drop = FALSE) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_fill_manual(
    values = c(
      "Stool / faeces" = "#66C2A5",
      "Blood / serum / urine" = "#FC8D62",
      "Tumour / mucosa" = "#8DA0CB",
      "Preclinical model" = "#E78AC3",
      "Multi-omics / literature" = "#A6D854"
    ),
    labels = specimen_labels,
    name = "Evidence source"
  ) +
  labs(
    title = "Evidence source composition",
    x = "Within-axis proportion",
    y = NULL
  ) +
  theme_postbio_fig2(base_size = 10) +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_blank()
  )

## Panel C
panelC_df <- dat_fig2 %>%
  count(axis, evidence_strength_factor, name = "n") %>%
  group_by(axis) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  mutate(axis = factor(axis, levels = rev(axis_levels)))

pC <- ggplot(panelC_df, aes(x = frac, y = axis, fill = evidence_strength_factor)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.72) +
  geom_text(
    aes(label = ifelse(frac >= 0.12, n, "")),
    position = position_stack(vjust = 0.5),
    size = 2.8,
    color = "black"
  ) +
  scale_y_discrete(labels = axis_labels_short, drop = FALSE) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_fill_manual(
    values = c(
      "1" = "#FFF7BC",
      "2" = "#FEC44F",
      "3" = "#FD8D3C",
      "4" = "#F03B20",
      "5" = "#BD0026"
    ),
    labels = strength_labels,
    name = "Evidence strength"
  ) +
  labs(
    title = "Evidence maturity profile",
    x = "Within-axis proportion",
    y = NULL
  ) +
  theme_postbio_fig2(base_size = 10) +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_blank()
  )

fig2 <- (pA / (pB | pC)) +
  plot_layout(
    heights = c(1.25, 1.0),
    widths = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 15, color = "black")
    )
  )

out_pdf <- file.path(out_dir, "Fig2_literature_landscape_v02.pdf")
out_png <- file.path(out_dir, "Fig2_literature_landscape_v02.png")

ggsave(
  filename = out_pdf,
  plot = fig2,
  width = 13.2,
  height = 9.4,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = out_png,
  plot = fig2,
  width = 13.2,
  height = 9.4,
  units = "in",
  dpi = 600,
  bg = "white"
)

write_csv(panelA_df, file.path(out_dir, "Fig2_panelA_bubble_input_v02.csv"))
write_csv(panelB_df, file.path(out_dir, "Fig2_panelB_specimen_composition_input_v02.csv"))
write_csv(panelC_df, file.path(out_dir, "Fig2_panelC_evidence_maturity_input_v02.csv"))

cat("[", Sys.time(), "] Figure 2 v02 saved:\n", sep = "")
cat("PDF: ", out_pdf, "\n", sep = "")
cat("PNG: ", out_png, "\n", sep = "")
cat("[", Sys.time(), "] Figure 2 v02 script finished.\n", sep = "")
