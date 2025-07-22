phip_convert <- function(
    data_long_path,
    sample_id = NULL,
    peptide_id = NULL,
    subject_id = NULL,
    timepoint = NULL,
    present = NULL,
    fold_change = NULL,
    counts_input = NULL,
    counts_hit = NULL,
    backend = NULL,
    peptide_library = TRUE
) {
  # ------------------------------------------------------------------
  # 1. db-backend: default to "duckdb" if user supplies nothing
  # ------------------------------------------------------------------
  backend_choices <- c("arrow", "duckdb", "memory")

  backend <- if (is.null(backend)) {
    "duckdb" # implicit default
  } else {
    match.arg(backend, choices = backend_choices)
  }

  # ------------------------------------------------------------------
  # 2. resolving the data_long_file path to absolute
  # ------------------------------------------------------------------
  ## check if the data_long_path provided
  .chk_cond(missing(data_long_path) || !nzchar(data_long_path),
            "'data_long_path' must be provided and non‑empty,
            no default is set.")

  ## check if this is a directory or a file
  info <- file.info(data_long_path)

  ## pre-check if exists
  .chk_cond(all(is.na(info)),
            sprintf("Path '%s' does not exist", data_long_path))

  if (!info$isdir || is.na(info$isdir)) {
    .chk_path(data_long_path,
              "data_long_path",
              extension = c("csv", "parquet", "parq", "pq"))
  }

  ## resolve the path to data_long_path
  cfg <- .resolve_paths(data_long_path = data_long_path,
                        backend = backend,
                        peptide_library = peptide_library)

  ## filter the NULLs
  cfg <- Filter(Negate(is.null), cfg)

  print(cfg)
  # ------------------------------------------------------------------
  # 3. mapping the column names
  # ------------------------------------------------------------------

  # ------------------------------------------------------------------
  # 4. create the phip_data object with different backends
  # ------------------------------------------------------------------
}

