## 11_build_Fig4_TME_context_v01.R
## Final Figure 4:
## Immune and microenvironmental context of postbiotic-related axes in TCGA CRC

rm(list = ls())

## =========================
## 0) Packages
## =========================
pkgs <- c(
  "readr", "dplyr", "ggplot2", "patchwork",
  "stringr", "forcats", "scales"
)

to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(lapply(pkgs, library, character.only = TRUE))

## =========================
## 1) Paths
## =========================
project_root <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

dir_input <- file.path(project_root, "results", "figure_build_inputs")
dir_fig   <- file.path(project_root, "results", "figures_main")

dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

path_fig4A <- file.path(dir_input, "fig4A_TME_heatmap_input_v01.csv")
path_fig4B <- file.path(dir_input, "fig4B_MSI_assoc_input_v01.csv")
path_fig4C <- file.path(dir_input, "fig4C_immunehot_assoc_input_v01.csv")

stopifnot(file.exists(path_fig4A))
stopifnot(file.exists(path_fig4B))
stopifnot(file.exists(path_fig4C))

## =========================
## 2) Read data
## =========================
df_A_raw <- readr::read_csv(path_fig4A, show_col_types = FALSE)
df_B_raw <- readr::read_csv(path_fig4B, show_col_types = FALSE)
df_C_raw <- readr::read_csv(path_fig4C, show_col_types = FALSE)

message("Loaded Figure 4 inputs:")
message("  A: ", nrow(df_A_raw), " rows")
message("  B: ", nrow(df_B_raw), " rows")
message("  C: ", nrow(df_C_raw), " rows")

## =========================
## 3) Helpers
## =========================

axis_map <- c(
  "SCFA_axis"         = "SCFA",
  "Trp_axis"          = "Tryptophan",
  "Polyamine_axis"    = "Polyamine",
  "Ferroptosis_axis"  = "Ferroptosis",
  "SCFA"              = "SCFA",
  "Tryptophan"        = "Tryptophan",
  "Polyamine"         = "Polyamine",
  "Ferroptosis"       = "Ferroptosis"
)

tme_map <- c(
  "CD8_T_cells"         = "CD8 T cells",
  "Tregs"               = "Tregs",
  "Th1_like"            = "Th1-like",
  "Th17_like"           = "Th17-like",
  "NK_cells"            = "NK cells",
  "M1_macrophages"      = "M1 macrophages",
  "M2_macrophages"      = "M2 macrophages",
  "Cytolytic_activity"  = "Cytolytic activity",
  "Stromal_like"        = "Stromal-like",
  "Stromal_score"       = "Stromal-like",
  "CD8 T cells"         = "CD8 T cells",
  "Tregs"               = "Tregs",
  "Th1-like"            = "Th1-like",
  "Th17-like"           = "Th17-like",
  "NK cells"            = "NK cells",
  "M1 macrophages"      = "M1 macrophages",
  "M2 macrophages"      = "M2 macrophages",
  "Cytolytic activity"  = "Cytolytic activity",
  "Stromal-like"        = "Stromal-like"
)

axis_order <- c("Tryptophan", "SCFA", "Ferroptosis", "Polyamine")
tme_order  <- c(
  "CD8 T cells", "NK cells", "Cytolytic activity",
  "Th1-like", "Tregs",
  "M1 macrophages", "M2 macrophages",
  "Th17-like", "Stromal-like"
)

axis_pal <- c(
  "Tryptophan" = "#7B6DCC",
  "SCFA"       = "#E76F51",
  "Ferroptosis"= "#2A9D8F",
  "Polyamine"  = "#4C78A8"
)

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop("None of these columns found: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

pretty_axis <- function(x) {
  out <- unname(axis_map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

pretty_tme <- function(x) {
  out <- unname(tme_map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p)        ~ "",
    p < 0.001       ~ "***",
    p < 0.01        ~ "**",
    p < 0.05        ~ "*",
    TRUE            ~ ""
  )
}

theme_postbio <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1, color = "black"),
      plot.subtitle = element_text(size = base_size - 1, color = "black"),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      plot.margin = margin(6, 8, 6, 8)
    )
}

## =========================
## 4) Prepare Panel A
## =========================
axis_col_A <- pick_col(df_A_raw, c("axis_pretty", "axis_label", "axis"))
tme_col_A  <- pick_col(df_A_raw, c("TME_signature_pretty", "TME_signature_label", "TME_signature"))
rho_col_A  <- pick_col(df_A_raw, c("rho", "cor", "correlation"))
fdr_col_A  <- pick_col(df_A_raw, c("FDR", "fdr", "padj", "adj_p"))

df_A <- df_A_raw %>%
  mutate(
    axis_pretty = pretty_axis(.data[[axis_col_A]]),
    tme_pretty  = pretty_tme(.data[[tme_col_A]]),
    rho         = as.numeric(.data[[rho_col_A]]),
    FDR         = as.numeric(.data[[fdr_col_A]]),
    sig         = sig_stars(FDR)
  ) %>%
  filter(axis_pretty %in% axis_order, tme_pretty %in% tme_order) %>%
  mutate(
    axis_pretty = factor(axis_pretty, levels = rev(axis_order)),
    tme_pretty  = factor(tme_pretty,  levels = tme_order)
  )

fill_lim <- max(abs(df_A$rho), na.rm = TRUE)
fill_lim <- max(fill_lim, 0.35)
fill_lim <- min(fill_lim, 0.85)

pA <- ggplot(df_A, aes(x = tme_pretty, y = axis_pretty, fill = rho)) +
  geom_tile(color = "#BDBDBD", linewidth = 0.6) +
  geom_text(aes(label = sig), size = 4.2, fontface = "bold", color = "black") +
  scale_fill_gradient2(
    low = "#4C78A8",
    mid = "white",
    high = "#D1495B",
    midpoint = 0,
    limits = c(-fill_lim, fill_lim),
    oob = squish,
    name = "Spearman\nrho"
  ) +
  labs(
    title = "Tumour microenvironment correlations",
    subtitle = "Asterisks indicate FDR-adjusted significance",
    x = NULL,
    y = NULL
  ) +
  theme_postbio(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 35, hjust = 1, vjust = 1
    ),
    legend.position = "right",
    aspect.ratio = 0.42
  )

## =========================
## 5) Prepare Panels B and C
## =========================
prep_diff_df <- function(df) {
  axis_col <- pick_col(df, c("axis_pretty", "axis_label", "axis"))
  diff_col <- pick_col(df, c("diff_median", "median_diff", "difference"))
  fdr_col  <- pick_col(df, c("FDR", "fdr", "padj", "adj_p"))
  
  df %>%
    mutate(
      axis_pretty = pretty_axis(.data[[axis_col]]),
      diff_median = as.numeric(.data[[diff_col]]),
      FDR         = as.numeric(.data[[fdr_col]]),
      sig         = sig_stars(FDR)
    ) %>%
    filter(axis_pretty %in% axis_order) %>%
    mutate(
      axis_pretty = factor(axis_pretty, levels = rev(axis_order))
    )
}

df_B <- prep_diff_df(df_B_raw)
df_C <- prep_diff_df(df_C_raw)

all_diff <- c(df_B$diff_median, df_C$diff_median)
xmin <- min(-0.05, min(all_diff, na.rm = TRUE) * 1.15)
xmax <- max(all_diff, na.rm = TRUE) * 1.22 + 0.03
x_limits <- c(xmin, xmax)

make_diff_plot <- function(df, title_text, xlab_text) {
  x_offset <- 0.03 * diff(range(x_limits))
  
  ggplot(df, aes(x = diff_median, y = axis_pretty, color = axis_pretty)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "#4D4D4D") +
    geom_segment(
      aes(x = 0, xend = diff_median, yend = axis_pretty),
      linewidth = 1.0,
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_point(size = 3.2, show.legend = FALSE) +
    geom_text(
      aes(
        x = diff_median + ifelse(diff_median >= 0, x_offset, -x_offset),
        label = sig
      ),
      color = "black",
      fontface = "bold",
      size = 4.2,
      show.legend = FALSE
    ) +
    scale_color_manual(values = axis_pal) +
    scale_x_continuous(
      limits = x_limits,
      expand = expansion(mult = c(0.01, 0.08))
    ) +
    labs(
      title = title_text,
      x = xlab_text,
      y = NULL
    ) +
    theme_postbio(base_size = 11) +
    theme(
      legend.position = "none"
    )
}

pB <- make_diff_plot(
  df = df_B,
  title_text = "MSI-H versus MSS/MSI-L",
  xlab_text = "Median difference in axis score\n(MSI-H minus MSS/MSI-L)"
)

pC <- make_diff_plot(
  df = df_C,
  title_text = "Immune-high versus immune-low",
  xlab_text = "Median difference in axis score\n(Immune-high minus immune-low)"
)

## =========================
## 6) Combine final figure
## =========================

add_main_title <- FALSE

fig4 <- (pA / (pB + pC)) +
  plot_layout(heights = c(1.25, 1.0)) +
  plot_annotation(
    title = if (add_main_title) {
      "Immune and microenvironmental context of postbiotic-related axes in TCGA colorectal cancer"
    } else {
      NULL
    },
    caption = "* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.caption = element_text(size = 9, hjust = 1, color = "black"),
      plot.tag = element_text(face = "bold", size = 15, color = "black")
    )
  )

## =========================
## 7) Save outputs
## =========================
out_base <- file.path(dir_fig, "Fig4_TME_context_TCGA_clean_v01")

ggsave(
  filename = paste0(out_base, ".png"),
  plot = fig4,
  width = 12.5,
  height = 9.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = paste0(out_base, ".pdf"),
  plot = fig4,
  width = 12.5,
  height = 9.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

message("Figure 4 saved to:")
message("  ", paste0(out_base, ".png"))
message("  ", paste0(out_base, ".pdf"))
