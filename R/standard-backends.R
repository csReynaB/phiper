#' Build a **memory-backed** `phip_data` object (no DuckDB / Arrow) with
#' standard worflow
#'
#' Reads up one long data frame directly into R and returns an in-memory
#' `phip_data`. If the `data_long_path` is a directory, it reads all files in
#' the directory, assuming they have the same file structure. Heavy validation
#' of matrix content is deliberately
#' **omitted**—that should be performed by the `phip_data` class validator.
#'
#' @param cfg       Named list from `resolve_paths()`; must contain one
#'   `data_long_path` element. Can be a directory or a file.
#' @param colmap    List of colnames mapped to the reserved `phip_data`
#'   colnames. Useful when the colnames in the file are named different than the
#'   colnames of the `phip_data` object.
#'
#' @return A `phip_data` object whose `data_long` slot is a tibble in RAM.
#' @keywords internal
.standard_read_memory_backend <- function(cfg, colmap) {
  ## checking if path is dir or file
  info <- file.info(cfg$data_long_path)

  # ---------- 1.  read the data files -----------------------------------------
  if (info$isdir) {
    # 1) list all files in the directory (parquet, csv)
    files <- list.files(
      path       = cfg$data_long_path,
      pattern    = "\\.(parquet|parq|pq|csv)$",
      full.names = TRUE
    )

    # 2) read each file with .auto_read() and row‑bind into one
    # data.frame/tibble
    out <- do.call(
      rbind,
      lapply(files, .auto_read)
    )
  } else {
    out <- .auto_read(cfg$data_long_path)
  }

  ## clear memory
  rm(mats)
  invisible(gc())

  # ---------- 1.  read the data files -----------------------------------------
  out <- .rename_to_standard(out, colmap)

  # ---------- 4.  finish ------------------------------------------------------
  tibble::as_tibble(out)
}

.rename_to_standard <- function(df, colname_map) {
  # invert: actual names --> standard names
  actual_to_std <- setNames(names(colname_map), unlist(colname_map))

  # only keep those that actually exist in df
  actual_to_std <- actual_to_std[names(actual_to_std) %in% colnames(df)]

  # rename
  names(df)[names(df) %in% names(actual_to_std)] <- actual_to_std[names(df)[names(df) %in% names(actual_to_std)]]

  df
}
