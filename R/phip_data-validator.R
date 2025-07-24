#' Internal validator for <phip_data> objects
#' @keywords internal
validate_phip_data <- function(x,
                               na_warn_thresh = 0.50, # warn if >50 % NA / zero
                               auto_expand = TRUE) { # fill missing grid?
  .check_pd(x) # existing helper
  .data <- rlang::.data # to silence the lintr and R CMD CHECK

  tbl <- x$data_long
  cols <- colnames(tbl)
  reserved <- c(
    "subject_id", "sample_id", "timepoint",
    "peptide_id", "present", "fold_change",
    "counts_input", "counts_hit"
  )

  ## ---------------------------------------------- 1  STRUCTURE
  is_long <- all(c("subject_id", "timepoint") %in% cols)
  need <- if (is_long) {
    c("subject_id", "timepoint", "sample_id", "peptide_id")
  } else {
    c("sample_id", "peptide_id")
  }

  miss <- setdiff(need, cols)
  .chk_cond(
    length(miss) > 0,
    sprintf(
      "Missing mandatory column(s): %s",
      paste(miss, collapse = ", ")
    )
  )

  # quick validate if at least two rows in the data
  n_rows <- tbl |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::pull(.data$n)

  .chk_cond(
    n_rows < 1,
    "The `data_long` has no rows! No peptides and/or subjects are specified."
  )

  ## ---------------------------------------------- 2  OUTCOME FAMILY
  have <- list(
    present     = "present" %in% cols,
    fold_change = "fold_change" %in% cols,
    raw_counts  = all(c("input_count", "hit_count") %in% cols)
  )
  k <- sum(unlist(have))
  .chk_cond(
    k == 0,
    "No outcome column found (need present /
            fold_change / raw_counts)."
  )

  ## --------------------------------------------- 3  RESERVED-NAME COLLISIONS
  overlap <- intersect(x$meta$extra_cols, reserved)

  .chk_cond(
    length(overlap) > 0,
    sprintf(
      "extra_cols overlap with reserved names: %s",
      paste(overlap, collapse = ", ")
    )
  )

  ## --------------------------------------------- 4  ATOMIC COLUMNS
  sample0 <- if (inherits(tbl, "tbl_lazy")) {
    # LIMIT 0 --> fetch only the schema, then drop Arrow classes
    tbl |>
      utils::head(0) |>
      dplyr::collect() |>
      as.data.frame()
  } else {
    as.data.frame(tbl)
  }

  bad_atomic <- names(sample0)[vapply(sample0, is.list, logical(1))]

  .chk_cond(
    length(bad_atomic) > 0,
    sprintf(
      "Non-atomic (list) columns found: %s",
      paste(bad_atomic, collapse = ", ")
    )
  )

  ## --------------------------------------------- 5  UNIQUENESS

  if (is_long) {
    dup <- tbl |>
      dplyr::count(.data$subject_id, .data$timepoint, .data$peptide_id) |>
      dplyr::filter(.data$n > 1) |>
      utils::head(1) |>
      dplyr::collect()

    .chk_cond(
      nrow(dup) > 0,
      "Each (subject_id, peptide_id, timepoint) must
              map to exactly one sample_id."
    )
  } else {
    dup <- tbl |>
      dplyr::count(.data$sample_id, .data$peptide_id) |>
      dplyr::filter(.data$n > 1) |>
      utils::head(1) |>
      dplyr::collect()
    .chk_cond(nrow(dup) > 0, "Duplicate sample_id values found.")
  }

  ## ---------------------------------------------- 6 VALUE RANGES
  if (have$present) {
    bad <- tbl |>
      dplyr::filter(!is.na(.data$present) & !(.data$present %in% c(0, 1))) |>
      utils::head(1) |>
      dplyr::collect()
    .chk_cond(nrow(bad) > 0, "`present` must be 0/1/NA.")
  }

  if (have$fold_change) {
    ok <- tbl |>
      # multiply by 1.0 to coerce integer --> double in the SQL
      dplyr::mutate(fold_change = .data$fold_change * 1.0) |>
      dplyr::summarise(all_finite = all(is.finite(.data$fold_change))) |>
      dplyr::pull(.data$all_finite)
    .chk_cond(!ok, "`fold_change` contains Inf/-Inf or NA.")
  }

  if (have$raw_counts) {
    neg <- tbl |>
      dplyr::filter(.data$input_count < 0 | .data$hit_count < 0) |>
      utils::head(1) |>
      dplyr::collect()
    .chk_cond(nrow(neg) > 0, "Raw counts must be non-negative.")
  }

  ## ---------------------------------------------- 7  SPARSITY WARNING
  if (have$present) {
    prop_na <- tbl |>
      dplyr::summarise(
        p = sum(dplyr::if_else(is.na(.data$present), 1, 0)) / dplyr::n()
      ) |>
      dplyr::pull(.data$p)

    if (prop_na > na_warn_thresh) {
      cli::cli_warn(sprintf("present column is %.0f %% NA.", prop_na * 100))
    }
  }

  if (have$raw_counts) {
    prop_zero <- tbl |>
      dplyr::summarise(
        p = sum(dplyr::if_else(.data$input_count == 0 &
          .data$hit_count == 0, 1, 0)) /
          dplyr::n()
      ) |>
      dplyr::pull(.data$p)

    if (prop_zero > na_warn_thresh) {
      cli::cli_warn(sprintf(
        "Both count columns are 0 in %.0f %% of rows.",
        prop_zero * 100
      ))
    }
  }

  ## ---------------------------------------------- 8  PEPTIDE-ID COVERAGE
  if (!is.null(x$peptide_library)) {
    missing_in_lib <- setdiff(
      tbl |> dplyr::distinct(.data$peptide_id) |> dplyr::pull(),
      x$peptide_library |> dplyr::distinct(.data$peptide_id) |> dplyr::pull()
    )
    missing_in_lib <- missing_in_lib[order(missing_in_lib)]

    .chk_cond(
      length(missing_in_lib) > 0,
      sprintf(
        "peptide_id not found in peptide_library (e.g. %s)",
        missing_in_lib[1]
      ),
      error = FALSE # emit warning instead of abort
    )
  }

  ## ---------------------------------------------- 9  COMPARISONS TABLE
  cmp <- x$comparisons
  if (nrow(cmp) > 0) {
    allowed_cmp_cols <- c("comparison", "group1", "group2", "variable")
    extra_cmp <- setdiff(colnames(cmp), allowed_cmp_cols)
    .chk_cond(
      length(extra_cmp) > 0,
      sprintf(
        "Unexpected columns in comparisons: %s",
        paste(extra_cmp, collapse = ", ")
      )
    )

    valid_labels <- if (is_long) {
      # which of the two exist?
      valid_cols <- dplyr::intersect(
        c("timepoint", "group"),
        colnames(tbl)
      )

      tbl |> # keep only the existing ones
        dplyr::select(dplyr::all_of(valid_cols)) |>
        dplyr::distinct() |> # unique row-wise combinations
        dplyr::collect() |> # pull the tiny result into R
        unlist(use.names = FALSE) |> # flatten to one vector
        unique() # final de-dup
    } else {
      tbl |>
        dplyr::distinct(.data$group) |>
        dplyr::pull(.data$group)
    }

    bad_cmp <- unique(c(
      setdiff(cmp$group1, valid_labels),
      setdiff(cmp$group2, valid_labels)
    ))
    .chk_cond(
      length(bad_cmp) > 0,
      sprintf(
        "comparisons refer to unknown group labels (e.g. %s)",
        bad_cmp[1]
      )
    )
  }

  ## ---------------------------------------------- 10 FULL GRID
  dims <- tbl |>
    dplyr::summarise(
      n_pep = dplyr::n_distinct(.data$peptide_id),
      n_smp = dplyr::n_distinct(.data$sample_id),
      n_obs = dplyr::n()
    ) |>
    dplyr::collect()

  expect <- dplyr::pull(dims, .data$n_pep) * dplyr::pull(dims, .data$n_smp)

  if (dplyr::pull(dims, .data$n_obs) != expect) {
    cli::cli_warn(sprintf(
      "Counts table was not a full peptide x sample
      grid (expanded %s to %s rows).",
      dplyr::pull(dims, .data$n_obs), expect
    ))
    if (isTRUE(auto_expand)) {
      # 1. create the full grid remotely (depends on the backend)
      make_full_grid <- function(tbl) {
        sample_tbl <- tbl |>
          dplyr::distinct(sample_id) |>
          dplyr::mutate(dummy = 1L)

        peptide_tbl <- tbl |>
          dplyr::distinct(peptide_id) |>
          dplyr::mutate(dummy = 1L)

        # CROSS JOIN via dummy -- works on any dplyr/dbplyr version
        sample_tbl |>
          dplyr::inner_join(peptide_tbl, by = "dummy") |>
          dplyr::select(-dummy) # clean-up
      }

      # 2. isolate the metadata that can be recycled
      non_recyclable_cols <- c(
        "present", "fold_change",
        "counts_control", "counts_hit"
      )

      recyclable_cols <- colnames(tbl)[colnames(tbl) %nin% non_recyclable_cols]
      non_recyclable_cols <- colnames(tbl)[colnames(tbl)
      %in% non_recyclable_cols]


      ## ---- 1. Define the key column -------------------------------------------
      # sample_id is the only identifier that must remain unique in tbl
      key_cols <- c("sample_id")

      ## ---- 2. Pick every other column as a candidate for metadata -------------
      # Anything that is NOT part of the key is a potential metadata column
      candidate_cols <- setdiff(colnames(tbl), key_cols)

      ## ---- 3. Test which candidate columns are constant within each sample ----
      # For every sample_id, check how many distinct values each column has.
      # If the count is always 1, that column is "recyclable".
      const_tbl <- tbl |>
        dplyr::group_by(sample_id) |>
        dplyr::summarise(
          dplyr::across(
            dplyr::all_of(candidate_cols),
            ~ dplyr::n_distinct(.x, na.rm = FALSE) == 1,   # TRUE if only one unique value
            .names = "const_{.col}"
          ),
          .groups = "drop"
        ) |>
        collect()

      # Recyclable columns  = columns TRUE for all sample_id rows
      recyclable_cols <- candidate_cols[
        sapply(
          candidate_cols,
          function(col) all(const_tbl[[paste0("const_", col)]])
        )
      ]

      ## define the inherently non-recycleable columns as special cases -->
      ## smart fallback for edge cases
      non_extra <- c("peptide_id", "present", "fold_change",
                     "counts_input", "counts_hit")

      recyclable_cols <- setdiff(recyclable_cols,
                                 non_extra)

      recyclable_cols <- c("sample_id", recyclable_cols)

      ## the only special column that is nor recyclable neither non-recyclable is
      ## the peptide_id --> it shouldn't be present in any of them

      # Non‑recyclable columns  = everything else
      non_recyclable_cols <- setdiff(candidate_cols, c(recyclable_cols, "peptide_id"))
      non_recyclable_cols <- intersect(colnames(tbl), non_recyclable_cols)

      sample_meta <- tbl |> # original table
        dplyr::select(dplyr::all_of(recyclable_cols)) |>
        dplyr::distinct() # one row per sample_id

      # 3. mark the original rows
      tbl <- tbl |> dplyr::mutate(.row_exists = 1L)

      # 3. join this in two steps
      tbl <- make_full_grid(tbl) |>
         dplyr::left_join(sample_meta, by = "sample_id") |>
        dplyr::left_join(
          dplyr::select(
            tbl, sample_id, peptide_id,
            non_recyclable_cols, .row_exists
          ),
          by = c("sample_id", "peptide_id")
        ) |>
        dplyr::mutate(across(
          dplyr::all_of(non_recyclable_cols),
          ~ dplyr::if_else(is.na(.row_exists), # introduced row?
            0, # …→ replace with 0
            .x
          ) # …→ keep original (NA or value)
        )) |>
        dplyr::select(-.row_exists) # clean-up


      # helper -----------------------------------------------------------------
      if(x$backend == "duckdb" || x$backend == "arrow") {


      .register_tbl <- function(tbl,
                                con,
                                name = "data_long",
                                materialise_table,
                                temporary  = TRUE) {

        if (materialise_table) {
          # ---- materialise with compute() -------------------------------------
          # * name        – quoted automatically
          # * temporary   – TRUE  → CREATE TEMPORARY TABLE  ...
          #                 FALSE → CREATE TABLE (persist in current schema)
          #
          # compute() returns a dplyr::tbl object lazily pointing to the new table

          out <- dplyr::compute(
            tbl,
            name       = name,
            temporary  = temporary,
            unique_indexes = NULL,   # pass your own if you need them
            analyze    = TRUE
          )

        } else {
          # ---- build / replace a VIEW ------------------------------------------
          qry_sql <- dbplyr::sql_render(tbl)

          DBI::dbExecute(
            con,
            sprintf(
              "CREATE OR REPLACE %sVIEW %s AS %s;",
              if (temporary) "TEMPORARY " else "",
              DBI::dbQuoteIdentifier(con, name),
              qry_sql
            )
          )
          out <- dplyr::tbl(con, name)
        }

        out
      }

      # after you produced `tbl` with the full grid logic
      x$data_long <- .register_tbl(
        tbl        = tbl,
        con        = x$meta$con,                 # open DuckDB connection
        name       = "data_long",
        materialise_table = x$meta$materialise_table    # "materialised" or "view"
      )

      } else {
        x$data_long <- tbl
      }
      x$meta$full_cross <- TRUE
    } else {
      cli::cli_warn(sprintf(
        "Counts table is not a full peptide x sample grid (%s / %s rows). ",
        dplyr::pull(dims, .data$n_obs), expect
      ))

      x$meta$full_cross <- FALSE
    }
  } else {
    x$meta$full_cross <- TRUE
  }

  invisible(x)
}
