#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(microDecon)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop(sprintf("Missing value after %s", flag))
  args[[idx + 1]]
}

has_flag <- function(flag) {
  flag %in% args
}

input_dir <- get_arg_value("--input_dir", "emu_results")
blank_name <- get_arg_value("--blank", "NFW")
output_prefix <- get_arg_value("--out_prefix", file.path(input_dir, "microdecon"))
renormalize <- has_flag("--renormalize")

if (!dir.exists(input_dir)) {
  stop(sprintf("Input directory not found: %s", input_dir))
}

pattern <- "\\.fastq_rel-abundance\\.tsv$"
files <- list.files(input_dir, pattern = pattern, full.names = TRUE)

if (length(files) == 0) {
  stop(sprintf("No files matching pattern '%s' found in %s", pattern, input_dir))
}

sample_name_from_path <- function(path) {
  sub("\\.fastq_rel-abundance\\.tsv$", "", basename(path))
}

make_taxa_string <- function(df) {
  # Build a semicolon-separated taxonomy string from any of the known EMU columns.
  # We keep it simple: concatenate available ranks (skipping empty/NA).
  rank_cols <- c(
    "superkingdom", "clade", "phylum", "class", "order", "family", "genus", "species",
    "subspecies", "species subgroup", "species group"
  )
  present <- intersect(rank_cols, colnames(df))
  if (length(present) == 0) {
    return(rep("", nrow(df)))
  }

  parts <- df[, present, drop = FALSE]
  parts <- lapply(parts, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    trimws(x)
  })
  parts <- as.data.frame(parts, check.names = FALSE, stringsAsFactors = FALSE)

  apply(parts, 1, function(row) {
    row <- row[row != ""]
    if (length(row) == 0) "" else paste(row, collapse = "; ")
  })
}

samples <- list()
# Accumulate taxonomy columns across files (for long-format output).
tax_meta <- NULL

standardize_missing <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

update_tax_meta <- function(tax_meta, df_tax) {
  df_tax$tax_id <- as.character(df_tax$tax_id)

  # Standardize missing values: NA means “unknown”.
  for (cn in setdiff(colnames(df_tax), "tax_id")) {
    df_tax[[cn]] <- standardize_missing(df_tax[[cn]])
  }

  if (is.null(tax_meta)) {
    return(df_tax)
  }

  # Add any new columns to tax_meta.
  new_cols <- setdiff(colnames(df_tax), colnames(tax_meta))
  for (cn in new_cols) {
    tax_meta[[cn]] <- NA_character_
  }

  # Add any missing columns to df_tax (so indexing works).
  missing_cols <- setdiff(colnames(tax_meta), colnames(df_tax))
  for (cn in missing_cols) {
    df_tax[[cn]] <- NA_character_
  }

  # Append new tax_ids.
  new_ids <- setdiff(df_tax$tax_id, tax_meta$tax_id)
  if (length(new_ids) > 0) {
    tax_meta <- rbind(
      tax_meta,
      df_tax[match(new_ids, df_tax$tax_id), colnames(tax_meta), drop = FALSE]
    )
  }

  # Fill NAs in tax_meta with info from df_tax.
  m <- match(df_tax$tax_id, tax_meta$tax_id)
  for (cn in setdiff(colnames(tax_meta), "tax_id")) {
    incoming <- df_tax[[cn]]
    if (all(is.na(incoming))) next
    cur <- tax_meta[[cn]]
    idx <- !is.na(incoming) & (is.na(cur[m]) | cur[m] == "")
    if (any(idx)) {
      cur[m[idx]] <- incoming[idx]
      tax_meta[[cn]] <- cur
    }
  }

  tax_meta
}

for (f in files) {
  sname <- sample_name_from_path(f)
  df <- read.delim(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  required_cols <- c("tax_id", "abundance")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("File %s is missing required columns: %s", f, paste(missing, collapse = ", ")))
  }

  ids <- as.character(df$tax_id)
  abund_raw <- suppressWarnings(as.numeric(df$abundance))
  if (any(is.na(abund_raw))) {
    bad <- sum(is.na(abund_raw))
    stop(sprintf("File %s has %d NA abundances after numeric conversion", f, bad))
  }

  # Aggregate in case EMU outputs duplicate tax_id rows.
  abund_t <- tapply(abund_raw, ids, sum, na.rm = TRUE)
  abund <- as.numeric(abund_t)
  names(abund) <- names(abund_t)
  samples[[sname]] <- abund

  tax_cols <- setdiff(colnames(df), c("tax_id", "abundance"))
  if (length(tax_cols) > 0) {
    df_tax <- df[, c("tax_id", tax_cols), drop = FALSE]
    df_tax$tax_id <- as.character(df_tax$tax_id)
    df_tax <- df_tax[!duplicated(df_tax$tax_id), , drop = FALSE]
    tax_meta <- update_tax_meta(tax_meta, df_tax)
  }
}

if (!(blank_name %in% names(samples))) {
  stop(sprintf(
    "Blank sample '%s' not found. Available samples: %s",
    blank_name,
    paste(sort(names(samples)), collapse = ", ")
  ))
}

all_ids <- sort(unique(unlist(lapply(samples, names), use.names = FALSE)))

abund_matrix <- sapply(samples, function(v) {
  out <- numeric(length(all_ids))
  names(out) <- all_ids
  idx <- match(names(v), all_ids)
  keep <- !is.na(idx)
  out[idx[keep]] <- v[keep]
  out
})

# Ensure matrix even if only one sample
abund_matrix <- as.matrix(abund_matrix)

sample_names <- colnames(abund_matrix)
other_samples <- setdiff(sample_names, blank_name)
other_samples <- sort(other_samples)
ordered_samples <- c(blank_name, other_samples)

if (length(other_samples) == 0) {
  stop("No non-blank samples found (only the blank file is present).")
}

if (is.null(tax_meta)) {
  tax_meta <- data.frame(tax_id = all_ids, check.names = FALSE, stringsAsFactors = FALSE)
}

# Ensure tax_meta covers all tax_ids.
missing_ids <- setdiff(all_ids, tax_meta$tax_id)
if (length(missing_ids) > 0) {
  pad <- data.frame(tax_id = missing_ids, check.names = FALSE, stringsAsFactors = FALSE)
  for (cn in setdiff(colnames(tax_meta), "tax_id")) {
    pad[[cn]] <- NA_character_
  }
  tax_meta <- rbind(tax_meta, pad[, colnames(tax_meta), drop = FALSE])
}

# Reorder taxonomy rows to match all_ids.
tax_meta <- tax_meta[match(all_ids, tax_meta$tax_id), , drop = FALSE]

# Keep a simple combined taxonomy string as well (helpful for microDecon "taxa" column).
taxa_vec <- make_taxa_string(
  data.frame(
    tax_meta[, setdiff(colnames(tax_meta), "tax_id"), drop = FALSE],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
)

input_df <- data.frame(
  OTU_ID = all_ids,
  abund_matrix[, ordered_samples, drop = FALSE],
  Taxa = taxa_vec,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

input_path <- paste0(output_prefix, "_input.tsv")
write.table(input_df, file = input_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# microDecon expects: OTU IDs, then blanks, then samples (grouped), then optional taxa.
# We treat all samples as one group, so numb.ind is a single value (n samples).
decontaminated <- decon(
  data = input_df,
  numb.blanks = 1,
  numb.ind = c(length(other_samples)),
  taxa = TRUE
)

extract_decon_table <- function(x) {
  if (is.data.frame(x)) return(x)

  if (is.list(x)) {
    # microDecon versions differ; try common component names first.
    candidates <- c(
      "final.table", "Final.table", "final_table", "Final_table",
      "decontaminated", "decontaminated.table", "decontaminated_table",
      "result", "results", "table", "data"
    )
    for (nm in candidates) {
      if (!is.null(x[[nm]]) && is.data.frame(x[[nm]])) return(x[[nm]])
    }

    # Fallback: first data.frame-like element.
    for (elt in x) {
      if (is.data.frame(elt)) return(elt)
    }
  }

  stop(
    "microDecon::decon() did not return a data.frame, and no data.frame element could be extracted. ",
    "Please inspect the object returned by decon() (e.g., with str()) to see where the table lives."
  )
}

decontaminated <- extract_decon_table(decontaminated)

# Some workflows prefer clamping tiny negatives to 0.
# Keep taxonomy columns intact; only adjust numeric sample columns.
num_cols <- setdiff(colnames(decontaminated), c("OTU_ID", "Taxa"))
for (cn in num_cols) {
  decontaminated[[cn]] <- pmax(0, as.numeric(decontaminated[[cn]]))
}

if (renormalize) {
  for (cn in num_cols) {
    s <- sum(decontaminated[[cn]], na.rm = TRUE)
    if (s > 0) decontaminated[[cn]] <- decontaminated[[cn]] / s
  }
}

output_path <- paste0(output_prefix, "_decontaminated.tsv")
desired_tax_order <- c(
  "species", "genus", "family", "order", "class", "phylum", "clade", "superkingdom",
  "subspecies", "species subgroup", "species group"
)

tax_cols_present <- intersect(desired_tax_order, colnames(tax_meta))
tax_cols_extra <- setdiff(setdiff(colnames(tax_meta), "tax_id"), desired_tax_order)
tax_cols_out_wide <- c(tax_cols_present, tax_cols_extra)

decon_wide <- decontaminated
colnames(decon_wide)[colnames(decon_wide) == "OTU_ID"] <- "tax_id"
decon_wide$tax_id <- as.character(decon_wide$tax_id)

# Drop the combined Taxa string from the wide output; keep explicit tax level columns.
if ("Taxa" %in% colnames(decon_wide)) {
  decon_wide$Taxa <- NULL
}

tax_out_wide <- tax_meta[, c("tax_id", tax_cols_out_wide), drop = FALSE]

decon_wide <- merge(
  decon_wide,
  tax_out_wide,
  by = "tax_id",
  all.x = TRUE,
  sort = FALSE
)

# Ensure no literal NA values are written for taxonomy columns.
for (cn in tax_cols_out_wide) {
  if (cn %in% colnames(decon_wide)) {
    decon_wide[[cn]] <- as.character(decon_wide[[cn]])
    decon_wide[[cn]][is.na(decon_wide[[cn]])] <- ""
  }
}

# Reorder columns: tax_id, blank+samples, then taxonomy levels.
sample_cols_wide <- setdiff(colnames(decon_wide), c("tax_id", tax_cols_out_wide))
decon_wide <- decon_wide[, c("tax_id", sample_cols_wide, tax_cols_out_wide), drop = FALSE]

write.table(decon_wide, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Also write a long-format table closer to the original EMU output:
# sample_id, tax_id, abundance, plus taxonomy columns.
sample_cols_out <- other_samples
long_list <- lapply(sample_cols_out, function(sid) {
  data.frame(
    sample_id = sid,
    tax_id = as.character(decontaminated$OTU_ID),
    abundance = as.numeric(decontaminated[[sid]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
long_df <- do.call(rbind, long_list)
long_df$row_id <- seq_len(nrow(long_df))

tax_out <- tax_meta
colnames(tax_out)[colnames(tax_out) == "tax_id"] <- "tax_id"
long_merged <- merge(long_df, tax_out, by = "tax_id", all.x = TRUE, sort = FALSE)
long_merged <- long_merged[order(long_merged$row_id), , drop = FALSE]
long_merged$row_id <- NULL

# Reorder columns to mimic EMU output style.
tax_cols_out <- setdiff(colnames(long_merged), c("sample_id", "tax_id", "abundance"))
long_merged <- long_merged[, c("sample_id", "tax_id", "abundance", tax_cols_out), drop = FALSE]

# Ensure there are no literal NA values in the long output.
for (cn in colnames(long_merged)) {
  if (is.numeric(long_merged[[cn]])) {
    long_merged[[cn]][is.na(long_merged[[cn]])] <- 0
  } else {
    long_merged[[cn]] <- as.character(long_merged[[cn]])
    long_merged[[cn]][is.na(long_merged[[cn]])] <- ""
  }
}

output_long_path <- paste0(output_prefix, "_decontaminated_long.tsv")
write.table(long_merged, file = output_long_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

message("Wrote microDecon input: ", input_path)
message("Wrote decontaminated table: ", output_path)
message("Wrote decontaminated long table: ", output_long_path)

if (renormalize) {
  message("Renormalization enabled (--renormalize): each sample column sums to 1.")
}
