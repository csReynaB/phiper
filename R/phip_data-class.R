################################################################################
## phip_data S3 class  ---------------------------------------------------------
################################################################################
#' @title Construct a **phip_data** object
#'
#' @description Creates a fully-validated S3 object that bundles the tidy
#'   PhIP-Seq counts (`data_long`), optional comparison definitions, a
#'   peptide-library annotation table, and backend metadata.  The function
#'   performs a minimal sanity check on *comparisons* and normalises the chosen
#'   storage
#' *backend* before returning the object (validation of the data itself
#'   happens via `validate_phip_data()` helper).
#'
#' @param data_long A tidy data frame (or `tbl_lazy`) with one row per
#'   `peptide_id` x `sample_id` combination. **Required.**
#' @param comparisons A data frame describing two-way contrasts
#'   (\code{comparison}, \code{group1}, \code{group2}, \code{variable});
#'   defaults to an empty tibble if \code{NULL}.
#' @param backend Character string specifying the storage engine to use:
#'   \code{"memory"}, \code{"duckdb"}, or \code{"arrow"}.  If \code{NULL}
#'   the implicit default is \code{"duckdb"}.
#' @param peptide_library A data frame with one row per \code{peptide_id}
#'   and its annotations.  If \code{NULL}, the package’s current default
#'   library is used.
#' @param meta Optional named list of metadata flags to pre-populate the
#'   \code{meta} slot (rarely needed by users).
#' @param auto_expand Logical. If `TRUE` and the input is **not** already the
#'   full Cartesian product of `sample_id` x `peptide_id`, the function fills in
#'   the missing combinations.
#'   * Columns that are constant within a `sample_id` (metadata) are duplicated
#'     to the newly created rows.
#'   * Measurement columns such as `fold_change`, `present`, raw counts, or any
#'     other non‑recyclable fields are initialised to 0.
#'   The expanded table replaces `data_long` in place.
#' @param materialise_table Logical (DuckDB and Arrow back‑ends only).
#'   If `FALSE` (default) the result is registered as a **view**.
#'   If `TRUE` the result is fully **materialised** and stored as a physical
#'   table, which speeds up repeated queries at the cost of extra memory/disk.

#' @return An object of class \code{"phip_data"}.
#'
#' @examples
#' \dontrun{
#' ## minimal constructor call
#' pd <- new_phip_data(
#'   data_long = tidy_counts,
#'   comparisons = NULL,
#'   backend = "duckdb",
#'   peptide_library = TRUE
#' )
#' }
#' ## list available backends
#' c("memory", "duckdb", "arrow")
#'
#' @export
new_phip_data <- function(data_long,
                          comparisons,
                          backend = c("memory", "duckdb", "arrow"),
                          peptide_library = TRUE,
                          auto_expand = TRUE,
                          materialise_table = TRUE,
                          meta = list()) {
  backend <- if (is.null(backend)) {
    "duckdb" # implicit default
  } else {
    match.arg(backend, choices = c("arrow", "duckdb", "memory"))
  }

  # quick sanity check
  if (!is.null(comparisons)) {
    .chk_cond(
      !inherits(comparisons, "data.frame"),
      "`comparisons` must be a data.frame (or NULL)."
    )
  } else {
    comparisons <- tibble::tibble( # empty placeholder
      comparison = character(),
      group1     = character(),
      group2     = character(),
      variable   = character()
    )
  }

  # --------------------------------------------------------------------------
  # Download the peptide metadata library
  # --------------------------------------------------------------------------
  if (peptide_library) {
    peptide_library <- get_peptide_meta()
  } else {
    peptide_library <- NULL
  }

  # --------------------------------------------------------------------------
  # Scan column names for automatic meta flags
  # --------------------------------------------------------------------------
  # colnames() works for tibble, tbl_dbi, and arrow_dplyr_query
  cols <- colnames(data_long)
  standard_cols <- c(
    "subject_id", "sample_id", "timepoint", "peptide_id",
    "present", "fold_change", "counts_input", "counts_hit"
  )

  meta$longitudinal <- all(c("timepoint", "sample_id") %in% cols)
  meta$exist <- "present" %in% cols
  meta$fold_change <- "fold_change" %in% cols
  meta$raw_counts <- all(c("counts_input", "counts_hit") %in% cols)
  meta$extra_cols <- cols[cols %nin% standard_cols]
  meta$peptide_con <- attr(peptide_library, "duckdb_con")
  meta$materialise_table <- materialise_table

  # check if the data is a full_cross --> that means all peptide x unique_id
  # combinations are present in the data
  is_full_cross <- function(tbl,
                            peptide = "peptide_id",
                            sample = "sample_id") {
    pep_sym <- rlang::sym(peptide)
    smp_sym <- rlang::sym(sample)

    dims <- tbl |>
      dplyr::summarise(
        n_pep = dplyr::n_distinct(!!pep_sym),
        n_smp = dplyr::n_distinct(!!smp_sym),
        n_obs = dplyr::n()
      ) |>
      dplyr::collect()

    dims$n_obs == dims$n_pep * dims$n_smp
  }

  # pick the correct column ----------------------------------------------
  meta$full_cross <- is_full_cross(data_long,
    peptide = "peptide_id",
    sample  = "sample_id"
  )

  # --------------------------------------------------------------------------

  # define the object
  obj <- structure(
    list(
      data_long       = data_long, # lazy tbl or tibble
      comparisons     = tibble::as_tibble(comparisons),
      backend         = backend,
      peptide_library = peptide_library,
      meta            = meta
    ),
    class = "phip_data"
  )

  # --------------------------------------------------------------------------
  # Validate the objects used to construct the phip_data
  # --------------------------------------------------------------------------
  # will stop on error, warn on non-fatal issues
  obj <- validate_phip_data(obj, auto_expand = auto_expand)

  # return clean validated phip_data
  return(obj)
}
