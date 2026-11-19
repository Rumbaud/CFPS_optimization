# ---------------------------------------------------------------
# Tecan Spark kinetic export -> tidy wide table  (CSV version)
# Multi-fluorophore version: handles files with multiple fluorophore
# sections (e.g. GFP_5ul, mCherry_5ul, GFP_2.5ul, mCherry_2.5ul).
#
# Sections sharing the same fluorophore prefix (e.g. all "GFP*") are
# merged into one output. Each fluorophore type has its own blank groups.
#
# Blank groups are now defined by WELL POSITION (A1, B2, C4, D9, ranges
# like "A1:A4" or "B2:D9"), not by SM sample labels.
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
# Helpers for well-based selection  (defined before CONFIG so wells() is usable)
# ============================================================================

# Parse a single well ID like "B2" -> list(row = 2, col = 2)
parse_well <- function(w) {
  m <- str_match(toupper(trimws(w)), "^([A-P])\\s*(\\d{1,2})$")
  if (is.na(m[1, 1])) stop("Bad well ID: '", w, "'")
  r <- match(m[1, 2], LETTERS)
  c <- as.integer(m[1, 3])
  if (c < 1 || c > 24) stop("Column out of range (1-24) in well: '", w, "'")
  list(row = r, col = c)
}

# Expand well specs into a character vector of well IDs.
# Accepts any mix of:
#   single wells            wells("B2", "C4", "D9")
#   a row stretch           wells("A1:A4")        -> A1 A2 A3 A4
#   a column stretch        wells("B3:F3")        -> B3 C3 D3 E3 F3
#   a rectangular block     wells("B2:D4")        -> B2 B3 B4 C2 ... D4
#   comma strings / vectors wells("A1, A2, B7"), wells(c("A1","A2"))
wells <- function(...) {
  specs <- unlist(list(...), use.names = FALSE)
  specs <- unlist(strsplit(specs, ","), use.names = FALSE)
  specs <- trimws(specs)
  specs <- specs[nzchar(specs)]
  
  out <- character(0)
  for (s in specs) {
    if (grepl(":", s, fixed = TRUE)) {
      parts <- strsplit(s, ":", fixed = TRUE)[[1]]
      if (length(parts) != 2) stop("Bad well range: '", s, "'")
      a <- parse_well(parts[1])
      b <- parse_well(parts[2])
      rr <- seq(min(a$row, b$row), max(a$row, b$row))
      cc <- seq(min(a$col, b$col), max(a$col, b$col))
      out <- c(out, paste0(rep(LETTERS[rr], each = length(cc)),
                           rep(cc, times = length(rr))))
    } else {
      p <- parse_well(s)
      out <- c(out, paste0(LETTERS[p$row], p$col))
    }
  }
  unique(out)
}

# Sort well IDs by row letter, then numeric column (A1, A2, ... B1, ...)
order_wells <- function(ws) {
  r <- match(substr(ws, 1, 1), LETTERS)
  c <- as.integer(substring(ws, 2))
  ws[order(r, c)]
}

# Shorten a long well list for console output
abbrev <- function(x, n = 8) {
  if (length(x) <= n) return(paste(x, collapse = ", "))
  paste0(paste(x[1:n], collapse = ", "), ", ... (", length(x), " total)")
}

# ============================================================================
# CONFIG  --  edit these between runs
# ============================================================================

# Full path to the raw Tecan CSV file
infile <- "data/raw/2026/261114_EK.csv"

# Parent directory where the dated output folder will be created.
output_root <- "data/processed/2026"

# Set to FALSE to write only the raw kinetics CSV and skip blank subtraction
# entirely (the blank_groups below are then ignored, but the fluorophore names
# are still needed — they are what finds the sections in the file).
do_normalize <- FALSE

# Blank groups per fluorophore, defined by WELL POSITION.
# Top-level names are the fluorophore prefixes (matched against section names
# like "GFP_5ul", "mCherry_2.5ul" — anything starting with the key matches).
#
# For each group:
#   blank   = wells whose per-timepoint mean is subtracted
#   samples = wells to subtract it from; use "rest" for every measured well
#             that is not a blank in any group of this fluorophore.
# Blank wells are dropped from the normalized output.
blank_groups_by_fluorophore <- list(
  GFP     = list(
    ),
  
  mCherry = list(
    )
)


# ============================================================================
# Helpers
# ============================================================================

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

# ---- 1. Plate layout (OPTIONAL: find the "<>" header row) ------------------
# If the Tecan method carries sample labels in the plate layout, they are used
# in the column names as "<label>__<well>". If the layout is absent or empty,
# columns are named by well ID alone - everything downstream works on wells.
layout <- data.frame(well = character(), sample = character(),
                     stringsAsFactors = FALSE)

header_row <- which(!is.na(raw[[1]]) & raw[[1]] == "<>")
if (length(header_row) == 0) {
  cat("No plate layout header ('<>') found - columns will be named by well.\n")
} else {
  header_row <- header_row[1]
  layout_rows <- (header_row + 1):(header_row + 16)
  layout_cols <- 2:25
  row_letters <- LETTERS[1:16]
  
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
  
  # Map wells -> sample labels and order replicates side-by-side.
  # With no usable layout, fall back to plate order and well-only column names.
  labelled <- layout %>%
    filter(well %in% names(well_data)) %>%
    mutate(
      sm_num  = as.integer(str_match(sample, "^SM(\\d+)_")[, 2]),
      rep_num = as.integer(str_match(sample, "_(\\d+)$")[, 2])
    ) %>%
    arrange(sm_num, rep_num)
  
  if (nrow(labelled) == 0) {
    cat("  No sample labels for these wells - naming columns by well ID.\n")
    labelled <- data.frame(well   = order_wells(names(well_data)),
                           sample = NA_character_,
                           stringsAsFactors = FALSE)
  }
  
  # Build numeric data frame
  out_num <- data.frame(Time_h = time_h, check.names = FALSE)
  for (k in seq_len(nrow(labelled))) {
    w   <- labelled$well[k]
    lab <- labelled$sample[k]
    col_name <- if (is.na(lab) || !nzchar(lab)) w else paste0(lab, "__", w)
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
  
  if (!do_normalize) {
    cat("  Normalization off - raw file only.\n")
    return(invisible(NULL))
  }
  
  # ---- Blank-subtracted output ---------------------------------------------
  norm_outfile <- file.path(out_dir,
                            paste0(file_stem, "_", fluoro_key,
                                   "_kinetics_normalized.csv"))
  
  blank_groups <- blank_groups_by_fluorophore[[fluoro_key]]
  sample_cols_all <- setdiff(colnames(out_num), "Time_h")
  
  # Column name -> well ID (column names are "<label>__<well>")
  col_well <- sub("^.*__", "", sample_cols_all)
  names(col_well) <- sample_cols_all
  present_wells <- unname(col_well)
  
  all_blank_wells <- unique(unlist(lapply(blank_groups, `[[`, "blank")))
  all_target_wells <- unlist(lapply(blank_groups, function(g) {
    if (is.character(g$samples) && length(g$samples) == 1 &&
        tolower(g$samples) == "rest") character(0) else g$samples
  }))
  all_target_wells <- unique(all_target_wells)
  
  missing_blanks  <- setdiff(all_blank_wells,  present_wells)
  missing_targets <- setdiff(all_target_wells, present_wells)
  if (length(missing_blanks) > 0) {
    warning(fluoro_key, ": Blank wells not found in data (skipped): ",
            abbrev(missing_blanks))
  }
  if (length(missing_targets) > 0) {
    warning(fluoro_key, ": Target wells not found in data (skipped): ",
            abbrev(missing_targets))
  }
  
  # Blank wells are dropped from the normalized table
  keep_cols <- sample_cols_all[!(col_well %in% all_blank_wells)]
  out_norm <- cbind(
    data.frame(Time_h = out_num$Time_h, check.names = FALSE),
    out_num[, keep_cols, drop = FALSE]
  )
  
  for (grp in blank_groups) {
    blank_w <- grp$blank
    targets <- grp$samples
    
    blank_cols <- sample_cols_all[col_well %in% blank_w]
    if (length(blank_cols) == 0) next
    
    blank_mean <- rowMeans(out_num[, blank_cols, drop = FALSE], na.rm = TRUE)
    
    if (is.character(targets) && length(targets) == 1 &&
        tolower(targets) == "rest") {
      target_cols <- keep_cols
      target_desc <- "all non-blank wells"
    } else {
      target_cols <- keep_cols[col_well[keep_cols] %in% targets]
      target_desc <- abbrev(unname(col_well[target_cols]))
    }
    
    cat(sprintf("  Blanks [%s]: averaging %d wells -> subtracting from %s\n",
                abbrev(unname(col_well[blank_cols])), length(blank_cols),
                target_desc))
    
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