## 31_make_Supplementary_Table_S4_gene_signature_membership.R
## Supplementary Table S4:
## Full gene membership of postbiotic-related axes and marker-based TME signatures
##
## Purpose:
##   1) Document all genes used for the four postbiotic-related host axes.
##   2) Document all genes used for the nine marker-based TME signatures.
##   3) Improve reproducibility for manuscript review.

rm(list = ls())

suppressPackageStartupMessages({
  pkgs <- c("dplyr", "readr", "stringr", "tibble")
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install) > 0) install.packages(to_install)
  invisible(lapply(pkgs, library, character.only = TRUE))
})

cat("[", Sys.time(), "] Supplementary Table S4 script started.\n", sep = "")

## =========================================================
## 1) Paths
## =========================================================

project_root <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

path_gene_sets <- file.path(
  project_root,
  "data", "processed", "annotations",
  "postbiotic_gene_sets_list_v01.rds"
)

path_tme_script <- file.path(
  project_root,
  "code", "R",
  "06_TCGA_TME_signatures.R"
)

dir_tables <- file.path(project_root, "results", "tables")
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)

out_csv <- file.path(
  dir_tables,
  "Supplementary_Table_S4_gene_signature_membership_v01.csv"
)

out_tsv <- file.path(
  dir_tables,
  "Supplementary_Table_S4_gene_signature_membership_v01.tsv"
)

out_docx <- file.path(
  dir_tables,
  "Supplementary_Table_S4_gene_signature_membership_v01.docx"
)

if (!file.exists(path_gene_sets)) {
  stop("Postbiotic gene-set RDS not found: ", path_gene_sets)
}

if (!file.exists(path_tme_script)) {
  stop("TME signature script not found: ", path_tme_script)
}

## =========================================================
## 2) Expected signatures
## =========================================================

expected_axes <- c(
  "SCFA_axis",
  "Trp_axis",
  "Polyamine_axis",
  "Ferroptosis_axis"
)

expected_tme <- c(
  "CD8_T_cells",
  "Tregs",
  "Th1_like",
  "Th17_like",
  "NK_cells",
  "M1_macrophages",
  "M2_macrophages",
  "Cytolytic_activity",
  "Stromal_like"
)

## =========================================================
## 3) Helper functions
## =========================================================

clean_axis_label <- function(x) {
  dplyr::recode(
    x,
    "SCFA_axis" = "SCFA",
    "Trp_axis" = "Tryptophan / AhR",
    "Polyamine_axis" = "Polyamine",
    "Ferroptosis_axis" = "Ferroptosis / redox",
    .default = x
  )
}

clean_tme_label <- function(x) {
  dplyr::recode(
    x,
    "CD8_T_cells" = "CD8 T cells",
    "Tregs" = "Regulatory T cells",
    "Th1_like" = "Th1-like",
    "Th17_like" = "Th17-like",
    "NK_cells" = "NK cells",
    "M1_macrophages" = "M1 macrophages",
    "M2_macrophages" = "M2 macrophages",
    "Cytolytic_activity" = "Cytolytic activity",
    "Stromal_like" = "Stromal-like",
    .default = stringr::str_replace_all(x, "_", " ")
  )
}

postbiotic_gene_role <- function(gene, axis_id) {
  dplyr::case_when(
    ## SCFA axis
    gene %in% c("FFAR2", "FFAR3", "HCAR2") ~
      "SCFA-sensing GPCR / receptor node",
    gene %in% c("SLC5A8", "SLC16A1", "SLC16A3") ~
      "Monocarboxylate / SCFA transport node",
    gene %in% c("HDAC1", "HDAC2", "HDAC3") ~
      "SCFA-sensitive epigenetic regulation node",
    gene %in% c("ACADS", "ACADSB", "ACAT1", "ECHS1", "HADH", "HADHB") ~
      "Short-chain fatty-acid and mitochondrial metabolic node",
    
    ## Tryptophan / AhR axis
    gene %in% c("AHR") ~
      "Aryl hydrocarbon receptor signalling node",
    gene %in% c("CYP1A1", "CYP1B1") ~
      "AhR-responsive xenobiotic / metabolic target",
    gene %in% c("IDO1", "TDO2") ~
      "Tryptophan catabolism enzyme",
    gene %in% c("KMO", "KYNU") ~
      "Kynurenine-pathway enzyme",
    gene %in% c("IL10", "IL22", "RORC") ~
      "Immune-regulatory / mucosal immunity node",
    
    ## Polyamine axis
    gene %in% c("ODC1") ~
      "Polyamine biosynthesis rate-limiting enzyme",
    gene %in% c("AMD1", "SRM", "SMS") ~
      "Polyamine biosynthesis enzyme",
    gene %in% c("SAT1", "PAOX", "SMOX") ~
      "Polyamine catabolism / oxidation enzyme",
    gene %in% c("MTAP") ~
      "Methionine salvage / polyamine-linked metabolic node",
    
    ## Ferroptosis / redox axis
    gene %in% c("GPX4", "SLC7A11") ~
      "Ferroptosis-suppressive antioxidant defence node",
    gene %in% c("ACSL4", "LPCAT3", "ALOX5", "ALOX12") ~
      "Lipid remodelling / lipid-peroxidation node",
    gene %in% c("FTH1", "FTL", "TFRC", "SLC40A1", "NCOA4") ~
      "Iron handling / ferritinophagy node",
    gene %in% c("NFE2L2", "KEAP1") ~
      "NRF2-KEAP1 redox-stress regulation node",
    gene %in% c("HMOX1") ~
      "Heme catabolism / iron-redox stress node",
    
    TRUE ~ paste0(
      "Curated host scoring node for ",
      clean_axis_label(axis_id),
      " axis"
    )
  )
}

tme_marker_category <- function(gene, sig_id) {
  paste0("Marker included in the ", clean_tme_label(sig_id), " signature")
}

## ---------------------------------------------------------
## Robust extraction of the TME list from 06_TCGA_TME_signatures.R
## ---------------------------------------------------------
## The previous version assumed object name = tme_signatures.
## This version searches for the list that contains CD8_T_cells = c(...),
## then walks backward to the nearest '<- list(' assignment and evaluates it.

extract_tme_signature_list <- function(script_path, expected_names) {
  
  lines <- readLines(script_path, warn = FALSE)
  
  ## Find a line inside the TME list.
  anchor_idx <- grep("\\bCD8_T_cells\\s*=\\s*c\\s*\\(", lines)
  
  if (length(anchor_idx) == 0) {
    stop(
      "Could not find 'CD8_T_cells = c(...)' in:\n",
      script_path,
      "\nSend the TME-signature definition block from 06_TCGA_TME_signatures.R."
    )
  }
  
  anchor_idx <- anchor_idx[1]
  
  ## Walk backward to the nearest object <- list( assignment.
  candidate_starts <- grep("<-\\s*list\\s*\\(", lines[seq_len(anchor_idx)])
  if (length(candidate_starts) == 0) {
    stop(
      "Found CD8_T_cells but could not find a preceding '<- list(' assignment in:\n",
      script_path
    )
  }
  
  start_idx <- candidate_starts[length(candidate_starts)]
  
  ## Extract object name.
  start_line <- lines[start_idx]
  object_name <- stringr::str_match(start_line, "^\\s*([A-Za-z0-9_.]+)\\s*<-\\s*list\\s*\\(")[, 2]
  
  if (is.na(object_name) || object_name == "") {
    stop("Could not parse object name from list assignment line:\n", start_line)
  }
  
  ## Find matching closing parenthesis for the list.
  balance <- 0
  end_idx <- NA_integer_
  
  for (i in start_idx:length(lines)) {
    line_i <- lines[i]
    
    ## This is robust enough for simple list(c(...)) definitions.
    open_n  <- stringr::str_count(line_i, fixed("("))
    close_n <- stringr::str_count(line_i, fixed(")"))
    
    balance <- balance + open_n - close_n
    
    if (i > start_idx && balance <= 0) {
      end_idx <- i
      break
    }
  }
  
  if (is.na(end_idx)) {
    stop("Could not identify the end of the TME signature list in: ", script_path)
  }
  
  list_code <- paste(lines[start_idx:end_idx], collapse = "\n")
  
  env <- new.env(parent = baseenv())
  
  tryCatch(
    eval(parse(text = list_code), envir = env),
    error = function(e) {
      stop(
        "Failed to parse/evaluate the extracted TME signature list.\n",
        "Object guessed: ", object_name, "\n",
        "Extracted code:\n", list_code, "\n\n",
        "Original error: ", e$message
      )
    }
  )
  
  out <- get(object_name, envir = env)
  
  if (!is.list(out) || is.null(names(out))) {
    stop("Extracted object '", object_name, "' is not a named list.")
  }
  
  missing <- setdiff(expected_names, names(out))
  if (length(missing) > 0) {
    stop(
      "Extracted TME list object '", object_name, "' but it is missing expected signatures: ",
      paste(missing, collapse = ", "),
      "\nExtracted names were: ",
      paste(names(out), collapse = ", ")
    )
  }
  
  out <- out[expected_names]
  out <- lapply(out, function(x) sort(unique(as.character(x))))
  
  attr(out, "object_name") <- object_name
  attr(out, "start_idx") <- start_idx
  attr(out, "end_idx") <- end_idx
  
  out
}

## =========================================================
## 4) Load postbiotic gene sets
## =========================================================

postbiotic_sets <- readRDS(path_gene_sets)

if (!is.list(postbiotic_sets) || is.null(names(postbiotic_sets))) {
  stop("postbiotic_gene_sets_list_v01.rds must be a named list.")
}

postbiotic_sets <- lapply(
  postbiotic_sets,
  function(x) sort(unique(as.character(x)))
)

missing_axes <- setdiff(expected_axes, names(postbiotic_sets))

if (length(missing_axes) > 0) {
  stop(
    "Missing expected postbiotic axes in RDS: ",
    paste(missing_axes, collapse = ", ")
  )
}

postbiotic_sets <- postbiotic_sets[expected_axes]

cat("[", Sys.time(), "] Loaded postbiotic axes:\n", sep = "")
print(sapply(postbiotic_sets, length))

## =========================================================
## 5) Extract TME signatures
## =========================================================

tme_sets <- extract_tme_signature_list(
  script_path = path_tme_script,
  expected_names = expected_tme
)

cat("[", Sys.time(), "] Extracted TME signatures from object: ",
    attr(tme_sets, "object_name"), "\n", sep = "")
cat("[", Sys.time(), "] TME signature block lines: ",
    attr(tme_sets, "start_idx"), " to ", attr(tme_sets, "end_idx"), "\n", sep = "")

print(sapply(tme_sets, length))

## =========================================================
## 6) Build Supplementary Table S4
## =========================================================

postbiotic_tbl <- tibble(
  signature_type = "Postbiotic-related host axis",
  signature_id = rep(names(postbiotic_sets), lengths(postbiotic_sets)),
  gene_symbol = unlist(postbiotic_sets, use.names = FALSE)
) %>%
  mutate(
    signature_label = clean_axis_label(signature_id),
    role_or_marker_category = postbiotic_gene_role(gene_symbol, signature_id),
    source_category = paste0(
      "Curated from KEGG, Reactome, MSigDB and CRC/postbiotic literature; ",
      "harmonised to HGNC symbols and fixed before downstream scoring"
    ),
    scoring_use = "Used to compute postbiotic-related axis scores",
    notes = paste0(
      "Expression z-standardized across samples; ",
      "axis score calculated as mean z-score across available genes"
    )
  )

tme_tbl <- tibble(
  signature_type = "Marker-based TME signature",
  signature_id = rep(names(tme_sets), lengths(tme_sets)),
  gene_symbol = unlist(tme_sets, use.names = FALSE)
) %>%
  mutate(
    signature_label = clean_tme_label(signature_id),
    role_or_marker_category = tme_marker_category(gene_symbol, signature_id),
    source_category = paste0(
      "Marker-based TME signature defined in 06_TCGA_TME_signatures.R; ",
      "used for TCGA immune/stromal correlation analysis"
    ),
    scoring_use = "Used to compute marker-based TME signature scores",
    notes = paste0(
      "Expression z-standardized across samples; ",
      "signature score calculated as mean z-score across available genes. ",
      "These signatures are marker summaries, not formal cell-deconvolution estimates."
    )
  )

stable_s4 <- bind_rows(postbiotic_tbl, tme_tbl) %>%
  mutate(
    signature_type = factor(
      signature_type,
      levels = c(
        "Postbiotic-related host axis",
        "Marker-based TME signature"
      )
    ),
    signature_id = factor(
      signature_id,
      levels = c(expected_axes, expected_tme)
    )
  ) %>%
  arrange(signature_type, signature_id, gene_symbol) %>%
  mutate(
    signature_type = as.character(signature_type),
    signature_id = as.character(signature_id)
  ) %>%
  select(
    `Signature type` = signature_type,
    `Signature ID` = signature_id,
    `Signature label` = signature_label,
    `Gene symbol` = gene_symbol,
    `Role / marker category` = role_or_marker_category,
    `Source category` = source_category,
    `Scoring use` = scoring_use,
    `Notes` = notes
  )

## =========================================================
## 7) QC checks
## =========================================================

if (any(is.na(stable_s4$`Gene symbol`) | stable_s4$`Gene symbol` == "")) {
  stop("Empty gene symbols found in Supplementary Table S4.")
}

dup_check <- stable_s4 %>%
  count(`Signature type`, `Signature ID`, `Gene symbol`) %>%
  filter(n > 1)

if (nrow(dup_check) > 0) {
  stop(
    "Duplicate gene memberships within the same signature were found:\n",
    paste(capture.output(print(dup_check)), collapse = "\n")
  )
}

summary_tbl <- stable_s4 %>%
  count(`Signature type`, `Signature ID`, `Signature label`, name = "n_genes") %>%
  arrange(`Signature type`, `Signature ID`)

cat("\nMembership summary:\n")
print(summary_tbl, n = Inf)

cat("\nTotal membership rows in Supplementary Table S4: ", nrow(stable_s4), "\n", sep = "")

## =========================================================
## 8) Write CSV and TSV
## =========================================================

write_csv(stable_s4, out_csv)
write_tsv(stable_s4, out_tsv)

cat("\n[WROTE] CSV: ", out_csv, "\n", sep = "")
cat("[WROTE] TSV: ", out_tsv, "\n", sep = "")

## =========================================================
## 9) Optional DOCX export
## =========================================================

if (requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("flextable", quietly = TRUE)) {
  
  doc <- officer::read_docx()
  
  doc <- officer::body_add_par(
    doc,
    "Supplementary Table S4. Full gene membership of postbiotic-related axes and marker-based tumour microenvironment signatures.",
    style = "heading 1"
  )
  
  doc <- officer::body_add_par(
    doc,
    paste0(
      "This table lists all genes used to compute the four postbiotic-related host transcriptomic axes ",
      "and the nine marker-based tumour microenvironment signatures. Postbiotic-axis genes were curated from ",
      "pathway databases and CRC/postbiotic literature, harmonised to HGNC symbols and fixed before downstream scoring. ",
      "TME signatures were treated as marker-based transcriptomic summaries rather than formal cell-deconvolution estimates."
    ),
    style = "Normal"
  )
  
  ft <- flextable::flextable(stable_s4)
  ft <- flextable::autofit(ft)
  ft <- flextable::fontsize(ft, size = 8, part = "all")
  ft <- flextable::fontsize(ft, size = 8.5, part = "header")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::align(ft, align = "left", part = "all")
  ft <- flextable::valign(ft, valign = "top", part = "all")
  
  doc <- flextable::body_add_flextable(doc, ft)
  
  doc <- officer::body_add_par(
    doc,
    paste0(
      "Abbreviations: CRC, colorectal cancer; HGNC, HUGO Gene Nomenclature Committee; ",
      "SCFA, short-chain fatty acid; TME, tumour microenvironment."
    ),
    style = "Normal"
  )
  
  print(doc, target = out_docx)
  
  cat("[WROTE] DOCX: ", out_docx, "\n", sep = "")
  
} else {
  cat(
    "\n[NOTE] officer/flextable not installed; DOCX was not written.\n",
    "CSV and TSV outputs were created successfully.\n",
    sep = ""
  )
}

cat("[", Sys.time(), "] Supplementary Table S4 script finished.\n", sep = "")
