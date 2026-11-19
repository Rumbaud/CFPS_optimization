# ---------------------------------------------------------------
# Tecan Spark kinetic export -> tidy wide table  (CSV version)
#
# Handles BOTH Tecan CSV layouts, auto-detected:
#
#   "block" layout  - one block per well down the sheet:
#                     well ID row, Time row, Temp row, Mean row, ...
#                     preceded by fluorophore section headers (e.g. "GFP").
#                     Sample labels optionally taken from the "<>" plate map.
#
#   "wide" layout   - a "Reading" header row with one column per well
#                     ("S0001 (B02)"), readings down the rows, elapsed
#                     time in the "avg. time [s]" column.
#
# Sections sharing the same fluorophore prefix are merged into one output.
# Each fluorophore type has its own blank groups, defined by WELL POSITION
# (A1, B2, C4, D9, ranges like "A1:A4" or "B2:D9").
#
# Input : <YYMMDD>.csv or <YYMMDD_N>.csv  (raw Tecan sheet exported as CSV)
# Output (one set of files per fluorophore type):
#   <stem>_<fluorophore>_kinetics_raw.csv
#   <stem>_<fluorophore>_kinetics_normalized.csv   (unless do_normalize = FALSE)
# ---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)

# ============================================================================
# Helpers for well-based selection  (defined before CONFIG so wells() is usable)
# ============================================================================

# Parse a single well ID like "B2" or "B02" -> list(row = 2, col = 2)
parse_well <- function(w) {
  m <- str_match(toupper(trimws(w)), "^([A-P])\\s*(\\d{1,2})$")
  if (is.na(m[1, 1])) stop("Bad well ID: '", w, "'")
  r <- match(m[1, 2], LETTERS)
  c <- as.integer(m[1, 3])
  if (c < 1 || c > 24) stop("Column out of range (1-24) in well: '", w, "'")
  list(row = r, col = c)
}

# Canonical well ID, unpadded: "B02" -> "B2"
norm_well <- function(w) {
  p <- parse_well(w)
  paste0(LETTERS[p$row], p$col)
}

# Expand well specs into a character vector of well IDs.
# Accepts any mix of:
#   single wells            wells("B2", "C4", "D9")
#   a row stretch           wells("A1:A4")        -> A1 A2 A3 A4
#   a column stretch        wells("B3:F3")        -> B3 C3 D3 E3 F3
#   a rectangular block     wells("B2:D4")        -> B2 B3 B4 C2 ... D4
#   comma strings / vectors wells("A1, A2, B7"), wells(c("A1","A2"))
# Zero-padded input ("B02") is accepted and normalised.
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
      out <- c(out, norm_well(s))
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
infile <- "data/raw/2026/261111_EK.csv"

# Parent directory where the dated output folder will be created.
output_root <- "data/processed/2026"

# Set to FALSE to write only the raw kinetics CSV and skip blank subtraction
# entirely (the blank_groups below are then ignored, but the fluorophore names
# are still needed - they are what finds the sections in the file).
do_normalize <- FALSE

# Blank groups per fluorophore, defined by WELL POSITION.
# Top-level names are the fluorophore prefixes, matched against the section
# header in a "block" file (e.g. "GFP_5ul") or, in a "wide" file, against the
# nearest text above the reading table (e.g. the method name "GFP kinetics RG").
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

# Classify a section name (e.g. "GFP_5ul", "GFP kinetics RG.skax") to a
# fluorophore key. Prefix match first, then substring. NA if nothing matches.
classify_fluorophore <- function(section_name, keys) {
  if (is.na(section_name)) return(NA_character_)
  for (k in keys) if (startsWith(section_name, k)) return(k)
  for (k in keys) if (str_detect(section_name, fixed(k))) return(k)
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

col1 <- raw[[1]]
col2 <- raw[[2]]
fluoro_keys <- names(blank_groups_by_fluorophore)

# ---- 1. Plate layout (OPTIONAL, "block" files only) ------------------------
# If the Tecan method carries sample labels in the "<>" plate map, they are
# used in the column names as "<label>__<well>". If the map is absent or
# empty, columns are named by well ID alone.
layout <- data.frame(well = character(), sample = character(),
                     stringsAsFactors = FALSE)

header_row <- which(!is.na(col1) & col1 == "<>")
if (length(header_row) > 0) {
  header_row <- header_row[1]
  layout_rows <- (header_row + 1):(header_row + 16)
  layout_cols <- 2:25
  for (i in seq_along(layout_rows)) {
    for (j in seq_along(layout_cols)) {
      val <- raw[layout_rows[i], layout_cols[j]]
      if (!is.na(val) && val != "None" && nzchar(val)) {
        layout <- rbind(layout,
                        data.frame(well   = paste0(LETTERS[i], j),
                                   sample = val,
                                   stringsAsFactors = FALSE))
      }
    }
  }
}
cat("Wells in plate map:", nrow(layout), "\n")

# ============================================================================
# 2. Parse sections  ->  list of list(name, time_s, well_data)
# ============================================================================

# Nearest text above `row` in column 1 that names a configured fluorophore.
nearest_fluoro_text <- function(row) {
  above <- seq_len(row - 1)
  cand <- above[!is.na(col1[above]) &
                  !is.na(vapply(col1[above],
                                function(s) classify_fluorophore(s, fluoro_keys),
                                character(1)))]
  if (length(cand) == 0) return(NA_character_)
  col1[cand[length(cand)]]
}

# ---- "wide" layout: a "Reading" header row, one column per well ------------
parse_wide_sections <- function() {
  hdr_rows <- which(!is.na(col1) & trimws(col1) == "Reading")
  sections <- list()
  
  for (hi in seq_along(hdr_rows)) {
    hr  <- hdr_rows[hi]
    end <- if (hi < length(hdr_rows)) hdr_rows[hi + 1] - 1 else nrow(raw)
    
    body    <- (hr + 1):end
    is_data <- !is.na(col1[body]) & grepl("^[0-9]+$", trimws(col1[body]))
    drows   <- body[is_data]
    if (length(drows) == 0) next
    
    hdr <- unlist(raw[hr, ], use.names = FALSE)
    well_pat  <- "\\(\\s*([A-P])\\s*(\\d{1,2})\\s*\\)"
    well_cols <- which(!is.na(hdr) & str_detect(hdr, well_pat))
    if (length(well_cols) == 0) next
    
    time_s <- as.numeric(raw[drows, 2])
    
    wd <- list()
    for (cc in well_cols) {
      m <- str_match(hdr[cc], well_pat)
      w <- paste0(m[1, 2], as.integer(m[1, 3]))
      wd[[w]] <- as.numeric(raw[drows, cc])
    }
    
    nm <- nearest_fluoro_text(hr)
    sections[[length(sections) + 1]] <-
      list(name = nm, time_s = time_s, well_data = wd, header_row = hr)
  }
  sections
}

# ---- "block" layout: one block per well ------------------------------------
parse_block_sections <- function() {
  well_pattern <- "^[A-P](?:[1-9]|1[0-9]|2[0-4])$"
  fluoro_section_pattern <-
    paste0("^(", paste(fluoro_keys, collapse = "|"), ")(_|$| )")
  section_rows <- which(!is.na(col1) &
                          str_detect(col1, fluoro_section_pattern) &
                          (is.na(col2) | !nzchar(col2)))
  if (length(section_rows) == 0) {
    stop("No fluorophore section headers found. Looked for rows starting with: ",
         paste(fluoro_keys, collapse = ", "))
  }
  
  block_starts <- which(!is.na(col1) & str_detect(col1, well_pattern))
  if (length(block_starts) == 0) stop("No per-well data blocks found.")
  
  block_section <- vapply(block_starts, function(br) {
    prior <- section_rows[section_rows < br]
    if (length(prior) == 0) NA_character_ else col1[prior[length(prior)]]
  }, character(1))
  
  sections <- list()
  for (sr in section_rows) {
    idx <- which(block_section == col1[sr])
    if (length(idx) == 0) next
    
    first_time <- as.numeric(unlist(raw[block_starts[idx[1]] + 1, -1]))
    time_s     <- first_time[!is.na(first_time)]
    
    wd <- list()
    for (i in idx) {
      br <- block_starts[i]
      w  <- norm_well(raw[br, 1])
      mr <- as.numeric(unlist(raw[br + 3, -1]))
      wd[[w]] <- mr[seq_along(time_s)]
    }
    sections[[length(sections) + 1]] <-
      list(name = col1[sr], time_s = time_s, well_data = wd, header_row = sr)
  }
  sections
}

file_layout <- if (any(!is.na(col1) & trimws(col1) == "Reading")) "wide" else "block"
cat("Detected file layout:", file_layout, "\n")

sections <- if (file_layout == "wide") parse_wide_sections() else parse_block_sections()
if (length(sections) == 0) stop("No readable data sections found in ", infile)

cat("Sections found:\n")
for (s in sections) {
  cat("  row", s$header_row, ":", ifelse(is.na(s$name), "<unnamed>", s$name),
      "-", length(s$well_data), "wells,", length(s$time_s), "readings\n")
}

unmatched <- vapply(sections,
                    function(s) is.na(classify_fluorophore(s$name, fluoro_keys)),
                    logical(1))
if (any(unmatched)) {
  warning("Section(s) with no matching fluorophore key - not processed. ",
          "Add a key to blank_groups_by_fluorophore matching the text above ",
          "the reading table.")
}

# ============================================================================
# 3. Process each fluorophore type
# ============================================================================
process_fluorophore <- function(fluoro_key) {
  cat("\n========== Processing fluorophore:", fluoro_key, "==========\n")
  
  mine <- Filter(function(s)
    identical(classify_fluorophore(s$name, fluoro_keys), fluoro_key), sections)
  if (length(mine) == 0) {
    cat("  No sections found for", fluoro_key, "- skipping.\n")
    return(invisible(NULL))
  }
  
  time_s <- mine[[1]]$time_s
  time_h <- time_s / 3600
  cat("  Sections:", length(mine), " Readings:", length(time_h), "\n")
  
  # Merge wells across sections of this fluorophore
  well_data <- list()
  for (s in mine) {
    for (w in names(s$well_data)) {
      v <- s$well_data[[w]]
      length(v) <- length(time_h)      # pad or trim to this fluorophore's time
      well_data[[w]] <- v
    }
  }
  cat("  Wells with data:", length(well_data), "\n")
  
  # Map wells -> sample labels where a plate map exists; otherwise plate order
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
  if (length(blank_groups) == 0) {
    cat("  No blank groups configured for", fluoro_key, "- raw file only.\n")
    return(invisible(NULL))
  }
  
  sample_cols_all <- setdiff(colnames(out_num), "Time_h")
  
  # Column name -> well ID (column names are "<label>__<well>" or "<well>")
  col_well <- sub("^.*__", "", sample_cols_all)
  names(col_well) <- sample_cols_all
  present_wells <- unname(col_well)
  
  all_blank_wells <- unique(unlist(lapply(blank_groups, `[[`, "blank")))
  all_target_wells <- unique(unlist(lapply(blank_groups, function(g) {
    if (is.character(g$samples) && length(g$samples) == 1 &&
        tolower(g$samples) == "rest") character(0) else g$samples
  })))
  
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