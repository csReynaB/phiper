#' @export
print.phip_data <- function(x, ...) {
  cat(cli::rule(
    left  = "<phip_data>",
    right = paste0("backend: ", x$backend)
  ), "\n\n")

  # ---- counts preview -------------------------------------------------------
  prev <- tryCatch(
    {
      tbl <- x$data_long
      if (inherits(tbl, c("tbl_dbi", "arrow_dplyr_query"))) {
        tbl |>
          utils::head(5) |>
          dplyr::collect()
      } else {
        utils::head(tbl, 5)
      }
    },
    error = function(e) tibble::tibble(.error = "<preview failed>")
  )

  cat(cli::col_cyan("counts (first 5 rows):"), "\n")
  print(prev)

  ## --- size summary ----------------------------------------------------------
  n_cols <- ncol(x$data_long)

  .data <- rlang::.data ## to silence the R CMD CHECK and lintr
  n_rows <- tryCatch(
    if (inherits(x$data_long, c("tbl_dbi", "arrow_dplyr_query"))) {
      x$data_long |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::pull(.data$n)
    } else {
      nrow(x$data_long)
    },
    error = function(e) NA_integer_
  )
  cat("\n")
  cat(sprintf(
    "table size: %s rows x %s columns\n\n",
    format(n_rows, big.mark = ","), n_cols
  ))

  # ---- contrasts ------------------------------------------------------------
  cat(cli::col_cyan("contrasts:"), "\n")
  cat(paste0(knitr::kable(x$comparisons, format = "simple"), collapse = "\n"))
  cat("\n\n")

  # ---- peptide library preview ---------------------------------------------
  cat(cli::col_cyan("peptide library preview (first 5 rows):"), "\n")

  lib <- x$peptide_library

  # If the DuckDB connection was closed, reopen the cached table
  if (inherits(lib, "tbl_dbi")) {
    if (!DBI::dbIsValid(lib$src$con)) {
      lib <- get_peptide_meta()
      x$peptide_library <- lib # keep the fresh handle for later prints
    }
  }

  show_cols <- intersect(c("peptide_id", "pos", "len_seq"), colnames(lib))

  lib_preview <- tryCatch(
    {
      lib |>
        utils::head(5) |> # LIMIT 5
        dplyr::collect() |> # into R
        dplyr::select(dplyr::all_of(show_cols))
    },
    error = function(e) tibble::tibble(.error = "<preview failed>")
  )

  print(lib_preview)

  # summary: extra-column count + overall dim ---------------------------------
  n_cols <- length(colnames(lib))
  extra_nc <- n_cols - length(show_cols)

  n_rows <- tryCatch(
    lib |> dplyr::summarise(n = dplyr::n()) |> dplyr::pull(.data$n),
    error = function(e) NA_integer_
  )

  cat(sprintf("... plus %d more columns\n\n", extra_nc))
  cat(sprintf(
    "library size: %s rows x %s columns\n\n",
    format(n_rows, big.mark = ","), n_cols
  ))


  # ---- meta flags -----------------------------------------------------------
  cat(cli::col_cyan("meta flags:"), "\n")
  for (nm in names(x$meta)) {
    val <- x$meta[[nm]]
    val_txt <- if (is.atomic(val)) {
      paste(val, collapse = ", ")
    } else {
      paste0("<", class(val)[1], ">")
    }
    cat(sprintf("  %-15s %s\n", paste0(nm, ":"), val_txt))
  }
  cat("\n")

  invisible(x)
}


################################################################################
## plain accessors (cause no S3 generics) --------------------------------------
################################################################################
.check_pd <- function(obj) {
  .chk_cond(
    !inherits(obj, "phip_data"),
    "`x` must be a <phip_data> object."
  )
}

#' @title Retrieve the main PhIP-Seq counts table
#'
#' @description Quick accessor for the `data_long` slot of a **phip_data**
#' object.
#'
#' @param x A valid `phip_data` object.
#'
#' @return A tibble or lazy table with one row per peptide × sample pair.
#' @export
get_counts <- function(x) {
  .check_pd(x)
  x$data_long
}

#' @title Retrieve the comparisons definition table
#'
#' @description Returns the two-way comparison specifications stored inside a
#' **phip_data** object.
#'
#' @inheritParams get_counts
#' @return A tibble with columns: `comparison`, `group1`, `group2`,
#'   and `variable`.
#' @export
get_comparisons <- function(x) {
  .check_pd(x)
  x$comparisons
}

#' @title Report the storage backend in use
#'
#' @description Tells you whether the counts are held in `"duckdb"`, `"arrow"`,
#' or `"memory"`.
#'
#' @inheritParams get_counts
#' @return A single string (`"duckdb"`, `"arrow"`, or `"memory"`).
#' @export
get_backend <- function(x) {
  .check_pd(x)
  x$backend
}

#' @title Retrieve the metadata list
#'
#' @description Accesses the `meta` slot, which holds flags such as whether the
#'  table is a full peptide × sample grid, the available outcome columns, etc.
#'
#' @inheritParams get_counts
#' @return A named list.
#' @export
get_meta <- function(x) {
  .check_pd(x)
  x$meta
}

#' @title Retrieve the peptide-library annotation table
#'
#' @description Quick accessor for the `peptide_library` slot of a **phip_data**
#' object.
#'
#' @inheritParams get_counts
#' @return A tibble with one row per `peptide_id` and associated annotation.
#' @export
get_peptide_library <- function(x) {
  .check_pd(x)
  x$peptide_library
}

################################################################################
## dplyr verb wrappers  --------------------------------------------------------
## (dplyr already provides the generics, we only add methods)
################################################################################

# helper to copy x, modify counts, and return
.modify_pd <- function(.data, new_counts) {
  .data$data_long <- new_counts
  .data
}

#' @importFrom dplyr filter
#' @exportS3Method filter phip_data
filter.phip_data <- function(.data, ..., .preserve = FALSE) {
  .modify_pd(
    .data,
    dplyr::filter(.data$data_long, ..., .preserve = .preserve)
  )
}

#' @importFrom dplyr select
#' @exportS3Method select phip_data
select.phip_data <- function(.data, ...) {
  .modify_pd(.data, dplyr::select(.data$data_long, ...))
}

#' @importFrom dplyr mutate
#' @exportS3Method mutate phip_data
mutate.phip_data <- function(.data, ...) {
  .modify_pd(.data, dplyr::mutate(.data$data_long, ...))
}

#' @importFrom dplyr arrange
#' @exportS3Method arrange phip_data
arrange.phip_data <- function(.data, ...) {
  .modify_pd(.data, dplyr::arrange(.data$data_long, ...))
}

#' @importFrom dplyr summarise
#' @exportS3Method summarise phip_data
summarise.phip_data <- function(.data, ..., .groups = NULL) {
  dplyr::summarise(.data$data_long, ..., .groups = .groups)
}

#' @importFrom dplyr collect
#' @exportS3Method collect phip_data
collect.phip_data <- function(x, ...) {
  dplyr::collect(x$data_long, ...)
}

################################################################################
## helper to close duckdb connection -------------------------------
################################################################################

#' @title Disconnect backend database connections
#'
#' @description Closes any open database connections held inside a
#' `phip_data` object-both the `data_long` backend (e.g. DuckDB,
#' Arrow) and the optional `peptide_library` connection.
#'
#' Call this when you want to release resources explicitly; it is
#' also safe to rely on garbage collection to close the connections
#' automatically at the end of the R session.
#'
#' @param x A valid `phip_data` object.
#'
#' @return The input `phip_data` object, invisibly, with its
#'   connections closed.
#'
#' @export

disconnect <- function(x) {
  # close data_long backend, if any
  if (!is.null(x$meta$con) && DBI::dbIsValid(x$meta$con)) {
    DBI::dbDisconnect(x$meta$con, shutdown = TRUE)
  }
  # close peptide library connection
  if (!is.null(x$meta$peptide_con) && DBI::dbIsValid(x$meta$peptide_con)) {
    DBI::dbDisconnect(x$meta$peptide_con, shutdown = TRUE)
  }

  invisible(x)
}
