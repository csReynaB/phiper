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
      "Counts table was not a full peptide × sample grid (%s → %s rows).",
      dplyr::pull(dims, .data$n_obs), expect
    ))

    if (isTRUE(auto_expand)) {
      tbl_expanded <- phip_expand_full_grid(
        tbl,
        key_col = "sample_id",
        id_col = "peptide_id",
        add_exist = TRUE, # <- new flag: keep binary column
        exist_col = "exist",
        fill_override = list(
          present      = 0L,
          fold_change  = 0,
          input_count  = 0L,
          hit_count    = 0L,
          counts_input = 0L,
          counts_hit   = 0L
        )
      )

      if (x$backend %in% c("duckdb", "arrow")) {
        x$data_long <- phip_register_tbl(
          tbl_expanded,
          con = x$meta$con,
          name = "data_long",
          materialise_table = x$meta$materialise_table
        )
      } else {
        x$data_long <- tbl_expanded
      }
      x$meta$full_cross <- TRUE
    } else {
      cli::cli_warn(sprintf(
        "Counts table is not a full peptide × sample grid (%s / %s rows).",
        dplyr::pull(dims, .data$n_obs), expect
      ))
      x$meta$full_cross <- FALSE
    }
  } else {
    x$meta$full_cross <- TRUE
  }

  invisible(x)
}

# 1) keep your existing table-based function but rename it internally:
.phip_expand_full_grid_impl <- function(tbl,
                                        key_col = "sample_id",
                                        id_col = "peptide_id",
                                        fill_override = NULL,
                                        add_exist = FALSE,
                                        exist_col = "exist") {
  stopifnot(all(c(key_col, id_col) %in% colnames(tbl)))
  .data <- rlang::.data
  key_sym <- rlang::sym(key_col)
  id_sym <- rlang::sym(id_col)
  row_exists_sym <- rlang::sym(".row_exists")

  # 0) tiny sample to detect types
  sample0 <- if (inherits(tbl, "tbl_lazy")) {
    tbl |>
      utils::head(0) |>
      dplyr::collect() |>
      as.data.frame()
  } else {
    as.data.frame(tbl)
  }

  # 1) robust cross (only key & id survive)
  make_cross <- function(tbl) {
    key_tbl <- tbl |>
      dplyr::distinct(!!key_sym) |>
      dplyr::transmute(!!key_sym, dummy = 1L)
    id_tbl <- tbl |>
      dplyr::distinct(!!id_sym) |>
      dplyr::transmute(!!id_sym, dummy = 1L)
    key_tbl |>
      dplyr::inner_join(id_tbl, by = "dummy") |>
      dplyr::select(-dummy) |>
      dplyr::select(!!key_sym, !!id_sym)
  }

  # 2) recyclable vs non-recyclable
  candidate_cols <- setdiff(colnames(tbl), c(key_col, id_col))
  recyclable_cols <- character(0)
  if (length(candidate_cols)) {
    const_tbl <- tbl |>
      dplyr::group_by(dplyr::across(dplyr::all_of(key_col))) |>
      dplyr::summarise(
        dplyr::across(dplyr::all_of(candidate_cols),
          ~ dplyr::n_distinct(.x, na.rm = FALSE) == 1L,
          .names = "const_{.col}"
        ),
        .groups = "drop"
      ) |>
      dplyr::collect()

    recyclable_cols <- candidate_cols[vapply(
      candidate_cols,
      function(col) isTRUE(all(const_tbl[[paste0("const_", col)]])),
      logical(1)
    )]
  }
  recyclable_cols <- setdiff(recyclable_cols, id_col)
  non_recyclable_cols <- setdiff(candidate_cols, recyclable_cols)

  sample_meta <- if (length(recyclable_cols)) {
    tbl |>
      dplyr::select(dplyr::all_of(c(key_col, recyclable_cols))) |>
      dplyr::distinct()
  } else {
    tbl |>
      dplyr::select(dplyr::all_of(key_col)) |>
      dplyr::distinct()
  }

  base_cells <- tbl |>
    dplyr::select(dplyr::all_of(c(key_col, id_col, non_recyclable_cols))) |>
    dplyr::mutate(!!row_exists_sym := 1L)

  cross_tbl <- make_cross(tbl)

  # 3) joins
  by_key <- if (utils::packageVersion("dplyr") >= "1.1.0") dplyr::join_by(!!key_sym) else setNames(key_col, key_col)
  by_both <- if (utils::packageVersion("dplyr") >= "1.1.0") dplyr::join_by(!!key_sym, !!id_sym) else setNames(c(key_col, id_col), c(key_col, id_col))

  out <- cross_tbl |>
    dplyr::left_join(sample_meta, by = by_key) |>
    dplyr::left_join(base_cells, by = by_both)

  # 4) fill introduced rows based on .row_exists
  common_measurement_names <- c(
    "present", "fold_change",
    "input_count", "hit_count",
    "counts_input", "counts_hit"
  )
  have_cols <- intersect(non_recyclable_cols, colnames(out))
  fill_candidates <- intersect(common_measurement_names, colnames(out))
  more_candidates <- setdiff(have_cols, fill_candidates)

  is_num <- function(x) inherits(x, c("numeric", "integer", "double", "integer64"))
  is_lgl <- function(x) is.logical(x)
  to_fill_num <- unique(c(
    intersect(fill_candidates, names(sample0)[vapply(sample0, is_num, logical(1))]),
    more_candidates[vapply(more_candidates, function(cn) is_num(sample0[[cn]]), logical(1))]
  ))
  to_fill_lgl <- unique(c(
    intersect(fill_candidates, names(sample0)[vapply(sample0, is_lgl, logical(1))]),
    more_candidates[vapply(more_candidates, function(cn) is_lgl(sample0[[cn]]), logical(1))]
  ))

  # custom overrides first
  if (!is.null(fill_override) && length(fill_override)) {
    for (nm in intersect(names(fill_override), colnames(out))) {
      out <- out |>
        dplyr::mutate(
          !!rlang::sym(nm) :=
            dplyr::if_else(is.na(!!row_exists_sym), fill_override[[nm]], .data[[nm]])
        )
      to_fill_num <- setdiff(to_fill_num, nm)
      to_fill_lgl <- setdiff(to_fill_lgl, nm)
    }
  }

  if (length(to_fill_num)) {
    out <- out |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(to_fill_num),
          ~ dplyr::if_else(is.na(!!row_exists_sym), 0, .x)
        )
      )
  }
  if (length(to_fill_lgl)) {
    out <- out |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(to_fill_lgl),
          ~ dplyr::if_else(is.na(!!row_exists_sym), FALSE, .x)
        )
      )
  }

  # 5) optional public "exist" column (0/1 integer), then drop sentinel
  if (isTRUE(add_exist)) {
    # avoid accidental overwrite
    exist_name <- exist_col
    if (exist_name %in% colnames(out)) {
      exist_name <- paste0(exist_name, "_grid")
      if (requireNamespace("cli", quietly = TRUE)) {
        cli::cli_warn(sprintf("Column '%s' exists; using '%s' for existence flag.", exist_col, exist_name))
      }
    }
    exist_sym <- rlang::sym(exist_name)
    out <- out |>
      dplyr::mutate(!!exist_sym := dplyr::if_else(is.na(!!row_exists_sym), 0L, 1L))
  }

  out |> dplyr::select(-!!row_exists_sym)
}

# 2) public generic
#' Expand to a full `sample_id × peptide_id` grid
#'
#' Creates the full Cartesian product of samples and peptides and joins
#' back per-sample metadata. For rows introduced by the expansion,
#' numeric/integer cell-level columns are filled with `0` (and logicals
#' with `FALSE`) unless overridden with `fill_override`.
#'
#' Works with local data frames and lazy `dbplyr` tables. The `phip_data`
#' method updates and returns the object (keeping laziness unless you
#' later `compute()` / `collect()` or use [phip_register_tbl()]).
#'
#' @section Method dispatch:
#' * `phip_expand_full_grid.data.frame()` – returns a (local) data frame/tibble.
#' * `phip_expand_full_grid.tbl_lazy()`   – returns a lazy table (still lazy).
#' * `phip_expand_full_grid.phip_data()`  – updates `x$data_long` and returns `x`.
#'
#' @param x Input object: a data frame / tibble, a lazy `dbplyr` table, or a
#'   `<phip_data>` object.
#' @param key_col Name of the sample identifier column. Default `"sample_id"`.
#' @param id_col  Name of the peptide identifier column. Default `"peptide_id"`.
#' @param fill_override Optional named list of fill values for **introduced** rows,
#'   e.g. `list(present = 0L, fold_change = NA_real_)`. Overrides the defaults.
#' @param ... Passed on to methods; currently unused.
#'
#' @return
#' * data.frame/tbl_lazy methods: object of the same general kind (local or lazy).
#' * phip_data method: the updated `<phip_data>` object (invisibly or visibly,
#'   depending on how you call it).
#'
#' @examples
#' \dontrun{
#' # local tibble
#' df_expanded <- phip_expand_full_grid(df_long)
#'
#' # lazy (DuckDB)
#' df_lazy <- dplyr::tbl(con, "data_long")
#' df_lazy2 <- phip_expand_full_grid(df_lazy)
#' dplyr::compute(df_lazy2, name = "data_long_full", temporary = TRUE)
#'
#' # <phip_data>: update in place & reassign
#' pd <- phip_expand_full_grid(pd, fill_override = list(fold_change = NA_real_))
#' }
#'
#' @seealso [phip_register_tbl()], [validate_phip_data()]
#' @export
phip_expand_full_grid <- function(x, ...) UseMethod("phip_expand_full_grid")

# 3) data.frame / tibble method
#' @rdname phip_expand_full_grid
#' @export
phip_expand_full_grid.data.frame <- function(x,
                                             key_col = "sample_id",
                                             id_col = "peptide_id",
                                             fill_override = NULL,
                                             add_exist = FALSE,
                                             exist_col = "exist",
                                             ...) {
  .phip_expand_full_grid_impl(x,
    key_col = key_col,
    id_col = id_col,
    fill_override = fill_override,
    add_exist = add_exist,
    exist_col = exist_col
  )
}

# 4) lazy DBI/dbplyr table method
#' @rdname phip_expand_full_grid
#' @export
phip_expand_full_grid.tbl_lazy <- function(x,
                                           key_col = "sample_id",
                                           id_col = "peptide_id",
                                           fill_override = NULL,
                                           add_exist = FALSE,
                                           exist_col = "exist",
                                           ...) {
  .phip_expand_full_grid_impl(x,
    key_col = key_col,
    id_col = id_col,
    fill_override = fill_override,
    add_exist = add_exist,
    exist_col = exist_col
  )
}

# 5) <phip_data> method — updates and returns the object
#' @rdname phip_expand_full_grid
#' @export
phip_expand_full_grid.phip_data <- function(x,
                                            key_col = "sample_id",
                                            id_col = "peptide_id",
                                            fill_override = NULL,
                                            add_exist = FALSE,
                                            exist_col = "exist",
                                            ...) {
  # expand (stays lazy)
  tbl_expanded <- .phip_expand_full_grid_impl(
    x$data_long,
    key_col       = key_col,
    id_col        = id_col,
    fill_override = fill_override,
    add_exist     = add_exist,
    exist_col     = exist_col
  )

  # compute grid completeness on the lazy object (small collect)
  dims <- tbl_expanded |>
    dplyr::summarise(
      n_pep = dplyr::n_distinct(.data[[id_col]]),
      n_smp = dplyr::n_distinct(.data[[key_col]]),
      n_obs = dplyr::n()
    ) |>
    dplyr::collect()
  expected <- dims$n_pep * dims$n_smp
  x$meta$full_cross <- (dims$n_obs == expected)

  # also store share of existing rows if we added the flag
  if (isTRUE(add_exist)) {
    exist_prop <- tbl_expanded |>
      dplyr::summarise(prop = mean(.data[[exist_col]] * 1.0)) |>
      dplyr::pull(.data$prop)
    x$meta$exist_prop <- as.numeric(exist_prop)
    x$meta$exist <- TRUE
  }

  # register back (keep your laziness/materialisation policy)
  if (x$backend %in% c("duckdb", "arrow")) {
    x$data_long <- phip_register_tbl(
      tbl_expanded,
      con = x$meta$con,
      name = "data_long",
      materialise_table = x$meta$materialise_table
    )
  } else {
    x$data_long <- tbl_expanded
  }

  x
}

#' Register a lazy table back to the DB as a TABLE or VIEW
#'
#' Convenience wrapper to either materialise a lazy pipeline via
#' `dplyr::compute()` (TABLE) or emit a `CREATE [TEMP] VIEW AS ...`
#' for engines like DuckDB. Returns a `dplyr::tbl()` pointing to the
#' created object.
#'
#' @param tbl A lazy table (e.g., from `dbplyr`).
#' @param con A `DBI` connection.
#' @param name Target name to create.
#' @param materialise_table If `TRUE`, create a TABLE via `compute()`. If `FALSE`,
#'   create a VIEW.
#' @param temporary If `TRUE`, create a TEMP table/view where supported.
#'
#' @return A lazy `dplyr::tbl` referencing the new TABLE/VIEW.
#' @examples
#' \dontrun{
#' lazy <- dplyr::tbl(con, "data_long") |> dplyr::filter(fold_change > 0)
#' phip_register_tbl(lazy, con, name = "data_long_pos", materialise_table = TRUE)
#' }
#' @export
phip_register_tbl <- function(tbl, con, name = "data_long",
                              materialise_table = TRUE, temporary = TRUE) {
  if (isTRUE(materialise_table)) {
    dplyr::compute(
      tbl,
      name = name,
      temporary = temporary,
      unique_indexes = NULL,
      analyze = TRUE
    )
  } else {
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
    dplyr::tbl(con, name)
  }
}
