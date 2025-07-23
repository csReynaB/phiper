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
    peptide_library = TRUE) {
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
  .chk_cond(
    missing(data_long_path) || !nzchar(data_long_path),
    "'data_long_path' must be provided and non‑empty,
            no default is set."
  )

  ## check if this is a directory or a file
  info <- file.info(data_long_path)

  ## pre-check if exists
  .chk_cond(
    all(is.na(info)),
    sprintf("Path '%s' does not exist", data_long_path)
  )

  if (!info$isdir || is.na(info$isdir)) {
    .chk_path(data_long_path,
      "data_long_path",
      extension = c("csv", "parquet", "parq", "pq")
    )
  }

  ## resolve the path to data_long_path
  cfg <- .resolve_paths(
    data_long_path = data_long_path,
    backend = backend,
    peptide_library = peptide_library
  )

  ## filter the NULLs
  cfg <- Filter(Negate(is.null), cfg)

  # ------------------------------------------------------------------
  # 3. mapping the column names
  # ------------------------------------------------------------------
  colname_map <- list(
    sample_id    = sample_id %||% "sample_id",
    peptide_id   = peptide_id %||% "peptide_id",
    subject_id   = subject_id %||% "subject_id",
    timepoint    = timepoint %||% "timepoint",
    present      = present %||% "present",
    fold_change  = fold_change %||% "fold_change",
    counts_input = counts_input %||% "counts_input",
    counts_hit   = counts_hit %||% "counts_hit"
  )

  # ------------------------------------------------------------------
  # 4. create the phip_data object with different backends
  # ------------------------------------------------------------------
  if (backend == "memory") {
    data_long <- .standard_read_memory_backend(cfg, colname_map)

    # returning the phip_data object
    new_phip_data(
      data_long = data_long,
      comparisons = NULL,
      peptide_library = cfg$peptide_library,
      backend = "memory"
    )
  } else if (backend == "duckdb") {
    # using the helper (a lot of code is repeated in duckdb and arrow, so i
    # decided to export it into a separate internale helper to reuse it)
    con <- .standard_read_duckdb_backend(cfg, colname_map)

    ## duckdb-specific code
    long <- dplyr::tbl(con, "data_long")

    # returning the phip_data object
    new_phip_data(
      data_long = long,
      comparisons = NULL,
      peptide_library = cfg$peptide_library,
      backend = "duckdb",
      meta = list(con = con)
    )
  } else if (backend == "arrow") {
    ## check dependency
    rlang::check_installed(c("arrow"), reason = "arrow backend")

    # same as up - use helper to create the data
    con <- .read_duckdb_backend(cfg, meta_list)

    # arrow-specific code, create tempdir to store the data
    arrow_dir <- withr::local_tempdir(
      pattern = sprintf(
        "phip_arrow_%s_",
        format(Sys.time(), "%Y%m%d%H%M%OS6")
      ),
      .local_envir = parent.frame()
    )

    # store the data as .parquet (more efficient than plain .csv)
    DBI::dbExecute(
      con,
      sprintf(
        "COPY final_long TO %s (FORMAT 'parquet', PER_THREAD_OUTPUT TRUE);",
        DBI::dbQuoteString(con, arrow_dir)
      )
    )

    # open as arrow dataset
    long <- arrow::open_dataset(arrow_dir)
    comps <- if (DBI::dbExistsTable(con, "comparisons")) {
      dplyr::tbl(con, "comparisons") |> dplyr::collect()
    }

    # returning the phip_data object
    new_phip_data(
      data_long = long,
      comparisons = comps,
      peptide_library = cfg$peptide_library,
      backend = "arrow",
      meta = list(parquet_dir = arrow_dir)
    )
  }
}
