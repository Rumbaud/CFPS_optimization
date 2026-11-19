# ---------------------------------------------------------------
# Tecan Spark kinetic export -> tidy wide table  (CSV version)
# Multi-fluorophore version: handles files with multiple fluorophore
# sections (e.g. GFP_5ul, mCherry_5ul, GFP_2.5ul, mCherry_2.5ul).
#
# Sections sharing the same fluorophore prefix (e.g. all "GFP*") are
# merged into one output. Each fluorophore type has its own blank groups.
#
# Input : <YYMMDD>.csv or <YYMMDD_N>.csv  (raw Tecan "Result sheet" exported as CSV)
# Output (one set of files per fluorophore type):
#   <stem>_<fluorophore>_kinetics_raw.csv
#   <stem>_<fluorophore>_kinetics_normalized.csv
#
# Where <stem> is the input filename without extension (e.g. "260414"),
# and <fluorophore> is e.g. "GFP" or "mCherry".
# ---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)

# ============================================================================
# CONFIG  --  edit these between runs
# ============================================================================

# Full path to the raw Tecan CSV file
infile <- "data/raw/2026/260818_EK.csv"

# Parent directory where the dated output folder will be created.
output_root <- "data/processed/2026"

# Helper: expand "SM2:SM4" → c("SM2","SM3","SM4")
sample_range <- function(x) {
  parts <- strsplit(x, ":")[[1]]
  prefix <- gsub("[0-9]+$", "", parts[1])
  from   <- as.integer(gsub("^\\D+", "", parts[1]))
  to     <- as.integer(gsub("^\\D+", "", parts[2]))
  paste0(prefix, seq(from, to))
}

# Blank groups per fluorophore.
# Top-level names are the fluorophore prefixes (matched against section names
# like "GFP_5ul", "mCherry_2.5ul" — anything starting with the key matches).
# Each value is a list of blank groups, same format as before.
blank_groups_by_fluorophore <- list(
  GFP     = list(
            list(blank   =   sample_range("SM4:SM6"),
                 samples = sample_range("SM1:SM3")),
            list(blank   =   sample_range("SM10:SM12"),
                 samples = sample_range("SM7:SM9")),
            list(blank   =   sample_range("SM16:SM18"),
                 samples = sample_range("SM13:SM15")),
            list(blank   =   sample_range("SM22:SM24"),
                 samples = sample_range("SM19:SM21")),
            list(blank   =   sample_range("SM28:SM30"),
                samples = sample_range("SM25:SM27")),
            list(blank   =   sample_range("SM34:SM36"),
                samples = sample_range("SM31:SM33")),
            list(blank   =   sample_range("SM40:SM42"),
                 samples = sample_range("SM37:SM39")),
            list(blank   =   sample_range("SM46:SM48"),
                 samples = sample_range("SM43:SM45")),
            list(blank   =   sample_range("SM52:SM54"),
                 samples = sample_range("SM49:SM51")))
)

# ============================================================================
# Helpers
# ============================================================================

sample_base <- function(x) sub("_(\\d+)(__.*)?$", "", x)

fmt_sci <- function(x) {
  ifelse(is.na(x), NA_character_, sprintf("%.2e", x))
}

# Classify a section name (e.g. "GFP_5ul") to a fluorophore key (e.g. "GFP").
# Returns NA if no key in blank_groups_by_fluorophore matches.
classify_fluorophore <- function(section_name, keys) {
  for (k in keys) {
    if (startsWith(section_name, k)) return(k)
  }
  NA_character_
}

# ============================================================================
# Output folder: use input filename stem (without extension) as folder name
# ============================================================================
file_stem <- tools::file_path_sans_ext(basename(infile))
out_dir   <- file.path(output_root, file_stem)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
cat("Output folder:", out_dir, "\n")

# ============================================================================
# Read the CSV as a raw character grid
# ============================================================================
raw <- read.csv(
  infile,
  header           = FALSE,
  colClasses       = "character",
  stringsAsFactors = FALSE,
  fill             = TRUE,
  check.names      = FALSE,
  na.strings       = c("", "NA")
)

# ---- 1. Plate layout (auto-detect: find the "<>" header row) ---------------
header_row <- which(!is.na(raw[[1]]) & raw[[1]] == "<>")
if (length(header_row) == 0) {
  stop("Could not find plate layout header row (looking for '<>' in column 1)")
}
header_row <- header_row[1]
layout_rows <- (header_row + 1):(header_row + 16)
layout_cols <- 2:25
row_letters <- LETTERS[1:16]

layout <- data.frame(well = character(), sample = character(),
                     stringsAsFactors = FALSE)
for (i in seq_along(layout_rows)) {
  for (j in seq_along(layout_cols)) {
    val <- raw[layout_rows[i], layout_cols[j]]
    if (!is.na(val) && val != "None" && nzchar(val)) {
      well <- paste0(row_letters[i], j)
      layout <- rbind(layout,
                      data.frame(well = well, sample = val,
                                 stringsAsFactors = FALSE))
    }
  }
}
cat("Wells in layout:", nrow(layout), "\n")

# ---- 2. Locate per-well kinetic blocks and tag with fluorophore ------------
# Per-well blocks start with a row whose first cell is a well ID (e.g. "I11"),
# followed by a Time row, Temp row, Mean row, StDev row, then per-read rows.
#
# Each block belongs to the most recent fluorophore section header above it.
# A section header is a row where col1 matches "<prefix>_..." (e.g. "GFP_5ul")
# and col2 is empty. We detect these by looking for rows that start with one
# of our configured fluorophore keys.
well_pattern <- "^[A-P](?:[1-9]|1[0-9]|2[0-4])$"
col1 <- raw[[1]]
col2 <- raw[[2]]

fluoro_keys <- names(blank_groups_by_fluorophore)
# Build a regex matching any section header: ^(GFP|mCherry)_...
fluoro_section_pattern <- paste0("^(", paste(fluoro_keys, collapse = "|"), ")(_|$)")
section_rows <- which(!is.na(col1) &
                        str_detect(col1, fluoro_section_pattern) &
                        (is.na(col2) | !nzchar(col2)))

if (length(section_rows) == 0) {
  stop("No fluorophore section headers found. Looked for rows starting with: ",
       paste(fluoro_keys, collapse = ", "))
}

cat("Fluorophore sections found:\n")
for (sr in section_rows) cat("  row", sr, ":", col1[sr], "\n")

block_starts <- which(!is.na(col1) & str_detect(col1, well_pattern))

# Tag each block with the section it belongs to (most recent section_row above)
block_section <- sapply(block_starts, function(br) {
  prior <- section_rows[section_rows < br]
  if (length(prior) == 0) return(NA_character_)
  col1[prior[length(prior)]]
})

# Time vector (assumed shared across all blocks — read from the first block)
first_time_row <- block_starts[1] + 1
time_vals_s <- as.numeric(unlist(raw[first_time_row, -1]))
time_vals_s <- time_vals_s[!is.na(time_vals_s)]
time_h <- time_vals_s / 3600
cat("Cycles (rows):", length(time_h), "\n")

# ---- 3. Process each fluorophore type --------------------------------------
process_fluorophore <- function(fluoro_key) {
  cat("\n========== Processing fluorophore:", fluoro_key, "==========\n")
  
  # Find blocks belonging to any section with this fluorophore key
  match_idx <- which(!is.na(block_section) &
                       vapply(block_section, function(s)
                         identical(classify_fluorophore(s, fluoro_keys),
                                   fluoro_key),
                         logical(1)))
  if (length(match_idx) == 0) {
    cat("  No blocks found for", fluoro_key, "- skipping.\n")
    return(invisible(NULL))
  }
  cat("  Blocks found:", length(match_idx), "\n")
  
  # Pull the Mean row from each block
  well_data <- list()
  for (idx in match_idx) {
    br <- block_starts[idx]
    well <- raw[br, 1]
    mean_row <- as.numeric(unlist(raw[br + 3, -1]))
    mean_row <- mean_row[seq_along(time_vals_s)]
    well_data[[well]] <- mean_row
  }
  cat("  Wells with data:", length(well_data), "\n")
  
  # Map wells -> sample labels and order replicates side-by-side
  labelled <- layout %>%
    filter(well %in% names(well_data)) %>%
    mutate(
      sm_num  = as.integer(str_match(sample, "^SM(\\d+)_")[, 2]),
      rep_num = as.integer(str_match(sample, "_(\\d+)$")[, 2])
    ) %>%
    arrange(sm_num, rep_num)
  
  if (nrow(labelled) == 0) {
    cat("  No wells in layout matched data wells - skipping.\n")
    return(invisible(NULL))
  }
  
  # Build numeric data frame
  out_num <- data.frame(Time_h = time_h, check.names = FALSE)
  for (k in seq_len(nrow(labelled))) {
    w   <- labelled$well[k]
    lab <- labelled$sample[k]
    col_name <- paste0(lab, "__", w)
    out_num[[col_name]] <- well_data[[w]]
  }
  cat("  Samples in output:", ncol(out_num) - 1, "\n")
  
  # ---- Write raw kinetics CSV ----------------------------------------------
  outfile <- file.path(out_dir,
                       paste0(file_stem, "_", fluoro_key, "_kinetics_raw.csv"))
  out_fmt <- data.frame(
    Time_h = sprintf("%.2f", out_num$Time_h),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  for (cn in setdiff(colnames(out_num), "Time_h")) {
    out_fmt[[cn]] <- fmt_sci(out_num[[cn]])
  }
  write.csv(out_fmt, outfile, row.names = FALSE)
  cat("  Wrote:", outfile, "\n")
  
  # ---- Blank-subtracted output ---------------------------------------------
  norm_outfile <- file.path(out_dir,
                            paste0(file_stem, "_", fluoro_key,
                                   "_kinetics_normalized.csv"))
  
  blank_groups <- blank_groups_by_fluorophore[[fluoro_key]]
  sample_cols_all <- setdiff(colnames(out_num), "Time_h")
  col_base <- sample_base(sample_cols_all)
  names(col_base) <- sample_cols_all
  
  all_blank_labels  <- unlist(lapply(blank_groups, `[[`, "blank"))
  all_target_labels <- unlist(lapply(blank_groups, `[[`, "samples"))
  present_labels    <- unique(col_base)
  
  missing_blanks  <- setdiff(all_blank_labels,  present_labels)
  missing_targets <- setdiff(all_target_labels, present_labels)
  if (length(missing_blanks) > 0) {
    warning(fluoro_key, ": Blank labels not found in data (skipped): ",
            paste(missing_blanks, collapse = ", "))
  }
  if (length(missing_targets) > 0) {
    warning(fluoro_key, ": Target sample labels not found in data (skipped): ",
            paste(missing_targets, collapse = ", "))
  }
  
  keep_cols <- sample_cols_all[!(col_base %in% all_blank_labels)]
  out_norm <- cbind(
    data.frame(Time_h = out_num$Time_h, check.names = FALSE),
    out_num[, keep_cols, drop = FALSE]
  )
  
  for (grp in blank_groups) {
    blank_lab <- grp$blank
    targets   <- grp$samples
    
    blank_cols <- sample_cols_all[col_base %in% blank_lab]
    if (length(blank_cols) == 0) next
    
    blank_mean <- rowMeans(out_num[, blank_cols, drop = FALSE], na.rm = TRUE)
    cat(sprintf("  Blanks [%s]: averaging %d reps -> subtracting from %s\n",
                paste(blank_lab, collapse = ", "), length(blank_cols),
                paste(targets, collapse = ", ")))
    
    target_cols <- keep_cols[col_base[keep_cols] %in% targets]
    if (length(target_cols) == 0) next
    
    for (tc in target_cols) {
      out_norm[[tc]] <- out_norm[[tc]] - blank_mean
    }
  }
  
  norm_fmt <- data.frame(
    Time_h = sprintf("%.2f", out_norm$Time_h),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  for (cn in setdiff(colnames(out_norm), "Time_h")) {
    norm_fmt[[cn]] <- fmt_sci(out_norm[[cn]])
  }
  write.csv(norm_fmt, norm_outfile, row.names = FALSE)
  cat("  Wrote:", norm_outfile, "\n")
}

# ---- 4. Run for every configured fluorophore -------------------------------
for (fk in fluoro_keys) {
  process_fluorophore(fk)
}

cat("\nDone.\n")