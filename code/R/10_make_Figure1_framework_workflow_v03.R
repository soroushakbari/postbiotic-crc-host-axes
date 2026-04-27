## 10_make_Figure1_framework_workflow_v03.R
## Figure 1: Conceptual framework and analytical workflow
## Polished stable ggplot2 version

rm(list = ls())

suppressPackageStartupMessages({
  pkgs <- c("ggplot2", "patchwork", "dplyr", "grid")
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install) > 0) install.packages(to_install)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

cat("[", Sys.time(), "] Figure 1 v03 script started.\n", sep = "")

project_root <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

dir_out <- file.path(project_root, "results", "figures_main")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

## =========================
## Colours
## =========================

col_scfa  <- "#1B8A78"
col_trp   <- "#5B54B7"
col_poly  <- "#D98222"
col_ferro <- "#D9573F"

col_scfa_light  <- "#E2F3F0"
col_trp_light   <- "#ECEAFB"
col_poly_light  <- "#FFF0DB"
col_ferro_light <- "#FDE8E2"

col_grey_dark  <- "#2F3A45"
col_grey_mid   <- "#6A737D"
col_grey_line  <- "#59636D"
col_grey_light <- "#F5F7F9"

axis_cols <- c(
  "SCFA" = col_scfa,
  "Tryptophan / AhR" = col_trp,
  "Polyamine" = col_poly,
  "Ferroptosis / redox" = col_ferro
)

theme_fig1 <- function(base_size = 10.5) {
  theme_void(base_size = base_size) +
    theme(
      text = element_text(color = "black", family = "sans"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 10, 8, 10)
    )
}

box_df <- function(x, y, label, fill, colour, id,
                   size = 3.0, fontface = "bold") {
  data.frame(
    x = x,
    y = y,
    label = label,
    fill = fill,
    colour = colour,
    id = id,
    size = size,
    fontface = fontface,
    stringsAsFactors = FALSE
  )
}

rect_df <- function(xmin, xmax, ymin, ymax, fill, colour, id) {
  data.frame(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    fill = fill,
    colour = colour,
    id = id,
    stringsAsFactors = FALSE
  )
}

label_layer <- function(df, pad = 0.34, radius = 0.18, linewidth = 0.45) {
  geom_label(
    data = df,
    aes(x = x, y = y, label = label, fill = fill, colour = colour),
    size = df$size,
    fontface = df$fontface,
    lineheight = 0.88,
    label.padding = unit(pad, "lines"),
    label.r = unit(radius, "lines"),
    linewidth = linewidth,
    show.legend = FALSE
  )
}

rect_layer <- function(df, linewidth = 0.45, alpha = 1) {
  geom_rect(
    data = df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = fill, colour = colour),
    linewidth = linewidth,
    alpha = alpha,
    show.legend = FALSE
  )
}

arrow_seg <- function(x, y, xend, yend, colour = col_grey_line,
                      linewidth = 0.62, linetype = "solid") {
  geom_segment(
    data = data.frame(x = x, y = y, xend = xend, yend = yend),
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    colour = colour,
    linewidth = linewidth,
    linetype = linetype,
    arrow = arrow(length = unit(0.12, "in"), type = "closed")
  )
}

arrow_curve <- function(x, y, xend, yend, colour = col_grey_line,
                        linewidth = 0.50, curvature = 0.20,
                        linetype = "solid") {
  geom_curve(
    data = data.frame(x = x, y = y, xend = xend, yend = yend),
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    colour = colour,
    linewidth = linewidth,
    curvature = curvature,
    linetype = linetype,
    arrow = arrow(length = unit(0.11, "in"), type = "closed")
  )
}

## =========================
## Panel A
## =========================

panelA_bg <- bind_rows(
  rect_df(0.25, 2.05, 3.15, 8.35, "#F8FAFC", "#D3DAE0", "bg_input"),
  rect_df(5.45, 7.55, 3.15, 8.35, "#FBFBFA", "#D3DAE0", "bg_host"),
  rect_df(8.55, 10.15, 3.15, 8.35, "#F8FAFC", "#D3DAE0", "bg_endpoint")
)

input_box <- box_df(
  1.15, 5.75,
  "Microbiome-linked inputs\n\nGut microbial communities\nDietary substrates\nPostbiotic-related metabolites",
  "white", "#71808A", "input",
  size = 3.05, fontface = "bold"
)

axis_boxes <- bind_rows(
  box_df(
    3.45, 7.75,
    "SCFA axis\nacetate, propionate, butyrate\nFFAR2/3, HCAR2, SLC5A8, HDAC",
    col_scfa_light, col_scfa, "axis_scfa", 2.82, "bold"
  ),
  box_df(
    3.45, 6.35,
    "Tryptophan / AhR axis\ntryptophan, kynurenine, indoles\nIDO1, TDO2, AHR",
    col_trp_light, col_trp, "axis_trp", 2.82, "bold"
  ),
  box_df(
    3.45, 4.95,
    "Polyamine axis\nputrescine, spermidine, spermine\nODC1, AMD1, SAT1, SMOX",
    col_poly_light, col_poly, "axis_poly", 2.82, "bold"
  ),
  box_df(
    3.45, 3.55,
    "Ferroptosis / redox axis\niron, GSH, lipid peroxides\nGPX4, SLC7A11, ACSL4, NFE2L2",
    col_ferro_light, col_ferro, "axis_ferro", 2.82, "bold"
  )
)

host_boxes <- bind_rows(
  box_df(
    6.50, 7.70,
    "Epithelial differentiation\nand barrier integrity",
    col_scfa_light, col_scfa, "host_1", 2.70, "plain"
  ),
  box_df(
    6.50, 6.72,
    "Immune activation\nor suppression",
    col_trp_light, col_trp, "host_2", 2.70, "plain"
  ),
  box_df(
    6.50, 5.74,
    "Proliferation and\nstress adaptation",
    col_poly_light, col_poly, "host_3", 2.70, "plain"
  ),
  box_df(
    6.50, 4.76,
    "Oxidative stress and\nferroptotic vulnerability",
    col_ferro_light, col_ferro, "host_4", 2.70, "plain"
  ),
  box_df(
    6.50, 3.78,
    "Stromal, vascular and\nmyeloid contexture",
    "#F1F3F5", "#68727B", "host_5", 2.70, "plain"
  )
)

endpoint_boxes <- bind_rows(
  box_df(9.35, 7.70, "Stage /\nprogression", "white", "#68727B", "end_1", 2.70, "plain"),
  box_df(9.35, 6.72, "Overall\nsurvival", "white", "#68727B", "end_2", 2.70, "plain"),
  box_df(9.35, 5.74, "MSI and immune\nphenotype", "white", "#68727B", "end_3", 2.60, "plain"),
  box_df(9.35, 4.76, "Therapy\nresponse", "white", "#68727B", "end_4", 2.70, "plain"),
  box_df(9.35, 3.78, "Biomarker\ndevelopment", "white", "#68727B", "end_5", 2.70, "plain")
)

microbe_points <- data.frame(
  x = c(0.65, 0.90, 1.20, 1.52, 0.78, 1.08, 1.42, 1.64, 0.70, 1.35),
  y = c(7.55, 7.95, 7.60, 7.85, 6.92, 6.62, 6.88, 6.45, 4.05, 3.75),
  colour = c(col_scfa, col_trp, col_poly, col_ferro,
             col_scfa, col_trp, col_poly, col_ferro,
             "#8DA0CB", "#66C2A5"),
  size = c(4.2, 3.7, 3.9, 4.1, 3.4, 3.5, 3.4, 3.7, 3.4, 3.2)
)

pA <- ggplot() +
  rect_layer(panelA_bg, linewidth = 0.55, alpha = 1) +
  geom_point(
    data = microbe_points,
    aes(x = x, y = y, colour = colour, size = size),
    alpha = 0.95,
    show.legend = FALSE
  ) +
  label_layer(input_box, pad = 0.34, linewidth = 0.36) +
  label_layer(axis_boxes, pad = 0.30, linewidth = 0.50) +
  label_layer(host_boxes, pad = 0.25, linewidth = 0.40) +
  label_layer(endpoint_boxes, pad = 0.24, linewidth = 0.36) +
  
  arrow_seg(2.00, 7.35, 2.48, 7.75, "#A7B0B7", 0.52) +
  arrow_seg(2.00, 6.42, 2.48, 6.35, "#A7B0B7", 0.52) +
  arrow_seg(2.00, 5.40, 2.48, 4.95, "#A7B0B7", 0.52) +
  arrow_seg(2.00, 4.35, 2.48, 3.55, "#A7B0B7", 0.52) +
  
  arrow_seg(4.62, 7.75, 5.42, 7.70, col_scfa, 0.78) +
  arrow_seg(4.62, 6.35, 5.42, 6.72, col_trp, 0.78) +
  arrow_seg(4.62, 4.95, 5.42, 5.74, col_poly, 0.78) +
  arrow_seg(4.62, 3.55, 5.42, 4.76, col_ferro, 0.78) +
  
  arrow_seg(7.52, 7.70, 8.55, 7.70, col_grey_line, 0.70) +
  arrow_seg(7.52, 6.72, 8.55, 6.72, col_grey_line, 0.70) +
  arrow_seg(7.52, 5.74, 8.55, 5.74, col_grey_line, 0.70) +
  arrow_seg(7.52, 4.76, 8.55, 4.76, col_grey_line, 0.70) +
  arrow_seg(7.52, 3.78, 8.55, 3.78, col_grey_line, 0.70) +
  
  arrow_curve(6.95, 7.42, 6.95, 6.05, "#A0A7AE", 0.42, curvature = -0.20, linetype = "dotted") +
  arrow_curve(6.02, 6.45, 6.02, 5.05, "#A0A7AE", 0.42, curvature = 0.20, linetype = "dotted") +
  arrow_curve(6.95, 5.45, 6.95, 4.10, "#A0A7AE", 0.42, curvature = -0.20, linetype = "dotted") +
  
  annotate(
    "text", x = 0.12, y = 9.25,
    label = "A", fontface = "bold", size = 7, hjust = 0
  ) +
  annotate(
    "text", x = 0.55, y = 9.25,
    label = "Conceptual framework of postbiotic-related host signalling in colorectal cancer",
    fontface = "bold", size = 4.6, hjust = 0
  ) +
  annotate(
    "text", x = 6.50, y = 8.72,
    label = "Host transcriptomic programmes\nin the colorectal tumour microenvironment",
    size = 3.35, lineheight = 0.94, fontface = "bold", colour = col_grey_dark
  ) +
  annotate(
    "text", x = 9.35, y = 8.72,
    label = "Clinical and\ntranslational endpoints",
    size = 3.30, lineheight = 0.94, fontface = "bold", colour = col_grey_dark
  ) +
  annotate(
    "text", x = 5.35, y = 2.58,
    label = "Effects are context-dependent and may be favourable or adverse depending on axis, compartment and disease state.",
    size = 3.10, fontface = "italic", colour = "#4C555E"
  ) +
  
  scale_fill_identity() +
  scale_colour_identity() +
  scale_size_identity() +
  coord_cartesian(xlim = c(0, 10.35), ylim = c(2.15, 9.45), clip = "off") +
  theme_fig1(base_size = 10.5)

## =========================
## Panel B
## =========================

workflow_boxes <- bind_rows(
  box_df(
    1.30, 6.25,
    "Focused evidence map\nTable 3 and Figure 2",
    col_scfa_light, col_scfa, "w1", 2.85, "bold"
  ),
  box_df(
    1.30, 4.95,
    "Host gene-set curation\n4 transcriptomic axes",
    col_trp_light, col_trp, "w2", 2.85, "bold"
  ),
  box_df(
    3.75, 5.60,
    "Transcriptomic scoring\nand external validation",
    "white", "#59636D", "w3", 3.05, "bold"
  ),
  box_df(
    3.25, 4.40,
    "TCGA CRC\nn = 614",
    col_scfa_light, col_scfa, "w4", 2.55, "bold"
  ),
  box_df(
    4.30, 4.40,
    "GSE39582\nn = 585",
    col_trp_light, col_trp, "w5", 2.55, "bold"
  ),
  box_df(
    6.10, 6.55,
    "Stage\nassociation",
    col_poly_light, col_poly, "w6", 2.55, "bold"
  ),
  box_df(
    6.10, 5.80,
    "Overall\nsurvival",
    "#EAF1FA", "#4C78A8", "w7", 2.55, "bold"
  ),
  box_df(
    6.10, 5.05,
    "TME\ncorrelation",
    col_scfa_light, col_scfa, "w8", 2.55, "bold"
  ),
  box_df(
    6.10, 4.30,
    "MSI-H vs\nMSS/MSI-L",
    col_trp_light, col_trp, "w9", 2.45, "bold"
  ),
  box_df(
    6.10, 3.55,
    "Immune-high vs\nimmune-low",
    col_ferro_light, col_ferro, "w10", 2.42, "bold"
  ),
  box_df(
    8.95, 5.05,
    "Translational interpretation\nand axis prioritisation",
    "#F6F8FA", "#59636D", "w11", 3.00, "bold"
  )
)

legend_df <- data.frame(
  x = c(1.10, 2.45, 4.10, 5.55),
  y = rep(2.55, 4),
  label = c("SCFA", "Tryptophan / AhR", "Polyamine", "Ferroptosis / redox"),
  colour = c(col_scfa, col_trp, col_poly, col_ferro),
  stringsAsFactors = FALSE
)

pB <- ggplot() +
  
  ## Axis-colour accent strip above the scoring box, not across the text
  geom_segment(aes(x = 3.05, xend = 4.45, y = 6.36, yend = 6.36),
               colour = col_scfa, linewidth = 1.45, lineend = "round") +
  geom_segment(aes(x = 3.05, xend = 4.45, y = 6.22, yend = 6.22),
               colour = col_trp, linewidth = 1.45, lineend = "round") +
  geom_segment(aes(x = 3.05, xend = 4.45, y = 6.08, yend = 6.08),
               colour = col_poly, linewidth = 1.45, lineend = "round") +
  geom_segment(aes(x = 3.05, xend = 4.45, y = 5.94, yend = 5.94),
               colour = col_ferro, linewidth = 1.45, lineend = "round") +
  
  label_layer(workflow_boxes, pad = 0.30, linewidth = 0.45) +
  
  arrow_seg(2.25, 6.25, 2.95, 5.78, col_grey_line, 0.62) +
  arrow_seg(2.25, 4.95, 2.95, 5.42, col_grey_line, 0.62) +
  
  arrow_seg(4.65, 5.60, 5.15, 5.60, col_grey_line, 0.62) +
  arrow_seg(5.15, 5.60, 5.38, 6.55, col_grey_line, 0.50) +
  arrow_seg(5.15, 5.60, 5.38, 5.80, col_grey_line, 0.50) +
  arrow_seg(5.15, 5.60, 5.38, 5.05, col_grey_line, 0.50) +
  arrow_seg(5.15, 5.60, 5.38, 4.30, col_grey_line, 0.50) +
  arrow_seg(5.15, 5.60, 5.38, 3.55, col_grey_line, 0.50) +
  
  arrow_seg(6.85, 5.05, 7.78, 5.05, col_grey_line, 0.65) +
  annotate(
    "text", x = 0.12, y = 7.28,
    label = "B", fontface = "bold", size = 7, hjust = 0
  ) +
  annotate(
    "text", x = 0.55, y = 7.28,
    label = "Study design of the present work",
    fontface = "bold", size = 4.6, hjust = 0
  ) +
  
  geom_point(
    data = legend_df,
    aes(x = x, y = y, colour = colour),
    size = 3.5,
    show.legend = FALSE
  ) +
  geom_text(
    data = legend_df,
    aes(x = x + 0.17, y = y, label = label),
    hjust = 0,
    size = 2.95,
    colour = col_grey_dark
  ) +
  
  scale_fill_identity() +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0, 10.15), ylim = c(2.30, 7.48), clip = "off") +
  theme_fig1(base_size = 10.5)

## =========================
## Compose and save
## =========================

fig1 <- pA / pB +
  plot_layout(heights = c(1.62, 1.00))

out_pdf <- file.path(dir_out, "Figure1_framework_workflow_v03.pdf")
out_png <- file.path(dir_out, "Figure1_framework_workflow_v03.png")

ggsave(
  filename = out_pdf,
  plot = fig1,
  width = 13.2,
  height = 12.0,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = out_png,
  plot = fig1,
  width = 13.2,
  height = 12.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(dir_out, "Figure1A_conceptual_framework_v03.pdf"),
  plot = pA,
  width = 13.2,
  height = 7.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(dir_out, "Figure1B_study_workflow_v03.pdf"),
  plot = pB,
  width = 13.2,
  height = 5.0,
  units = "in",
  device = cairo_pdf
)

cat("[", Sys.time(), "] Figure 1 v03 saved:\n", sep = "")
cat("PDF: ", out_pdf, "\n", sep = "")
cat("PNG: ", out_png, "\n", sep = "")
cat("[", Sys.time(), "] Figure 1 v03 script finished.\n", sep = "")
