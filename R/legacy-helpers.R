# ------------------------------------------------------------------------------
#  helper: fastest possible CSV reader with delimiter sniffing -
#   variation on the Carlos's function --> added .parquet support
# ------------------------------------------------------------------------------
.auto_read <- function(path,
                       ...) {

  base <- basename(path)

  ext <- strsplit(base, "\\.", fixed = FALSE)[[1]][-1]
  ext <- paste0(tolower(ext), collapse = ".")

  ## ------------------------------------------------------------------ ##
  ##                1.  PARQUET branch                                  ##
  ## ------------------------------------------------------------------ ##
  if (ext %in% c("parquet", "parq", "pq")) {
    ## to avoid additional dependencies load the data using arrow, which
    ## which is already listed in the dependencies
    rlang::check_installed("arrow")

    arrow::read_parquet(path,
                        as_data_frame = TRUE,
                        ...)
  } else {

    ## ------------------------------------------------------------------ ##
    ##                2.  CSV / TSV branch (original)                     ##
    ## ------------------------------------------------------------------ ##
    hdr <- readLines(path, n = 1L, warn = FALSE)

    n_comma <- lengths(regmatches(hdr, gregexpr(",", hdr, fixed = TRUE)))
    n_semi  <- lengths(regmatches(hdr, gregexpr(";", hdr, fixed = TRUE)))
    sep     <- if (n_semi > n_comma) ";" else ","

    if (requireNamespace("data.table", quietly = TRUE)) {
      rlang::check_installed("data.table")
      data.table::fread(path,
                        sep          = sep,
                        data.table   = FALSE,
                        check.names  = FALSE,
                        showProgress = FALSE,
                        ...)
    } else {
      utils::read.csv(path,
                      header          = TRUE,
                      sep             = sep,
                      check.names     = FALSE,
                      stringsAsFactors = FALSE,
                      ...)
    }
  }
}

#' @title Prepare sample metadata and (optional) comparisons for legacy import
#'
#' @description Reads the `samples_file`, optionally reads `comparisons_file`
#'   and `timepoints_file`, verifies that any dummy-coded grouping columns are
#'   mutually exclusive (row sum == 1), and collapses them into a single `group`
#'   column. When `timepoints_file` present, it merges tha `samples` with
#'   `timepoints` and adds an person identifier and timepoint variable.
#'
#' @param samples_file      Absolute path to the samples CSV/Parquet.
#' @param comparisons_file  Absolute path to comparisons CSV/Parquet, or `NULL`.
#' @param timepoints_file   Absolute path to timepoints CSV/Parquet, or `NULL`.
#' @param extra_cols        Character vector of extra metadata columns to keep.
#'
#' @return A list with elements `samples`, `comparisons` (or `NULL`), and
#'   `extra_cols` (possibly augmented).
#' @keywords internal
.legacy_prepare_metadata <- function(samples_file,
                                    comparisons_file = NULL,
                                    timepoints_file  = NULL,
                                    extra_cols = character()) {

  # ---- samples -------------------------------------------------------------
  samples <- .auto_read(samples_file)        # small table
  names(samples)[1] <- "sample_id"           # rename first var

  # this is my personal discussion with myself, maybe somebody will find this
  # entertaining in the future:

  # ----
  # the legacy workflow has an important limitation: the comparisons specified
  # in the comparisons_file are suitable only for the longitudinal data -->
  # this means that both timepoints_file and comparisons_file have to be
  # provided at the same time --> log error if only one provided

  ## ==> this is actually wrong, solution: see long comment around line 180 with
  ## XXX at the beginning

  # XXX
  # soooo, if i get it right, the groups in Carlos's script have to be defined
  # in the comparisons file and in the metadata as columns; so my idea to
  # solve this and allow comparisons for both long and cross-sectional data is
  # to first pull all unique values from the group1 and group2 cols in the
  # comparisons file, then select respective column from the metadata table.
  # They are dummy-coded so we can merge them into one single column: group
  # ----

  # ---- comparisons (optional) ---------------------------------------------
  if (is.null(comparisons_file)) {
    comparisons <- NULL
  } else {
    comparisons <- .auto_read(comparisons_file)

    # collect dummy-coded column names referenced in comparisons
    dummy_cols <- unique(c(comparisons$group1, comparisons$group2))

    # sanity-check: each row must belong to exactly one dummy group
    .chk_cond(
      any(rowSums(samples[, dummy_cols]) != 1),
      "Grouping columns in samples_file must be
      mutually exclusive (row sum != 1)."
    )

    # collapse dummies --> ‘group’
    which_max <- max.col(samples[, dummy_cols], ties.method = "first")
    samples$group <- names(samples[, dummy_cols])[which_max]
    samples       <- samples[, !(names(samples) %in% dummy_cols)]

    comparisons$variable <- "group"
  }

  ## delete the columns: only sample_id and group are allowed to stay
  keep <- colnames(samples) %in% c("sample_id", "group", extra_cols)
  samples <- samples[keep]

  # ---- time-points (optional) ---------------------------------------------
  if (is.null(timepoints_file)) {
    timepoints <- NULL
  } else {
    tp_wide <- .auto_read(timepoints_file)

    tp_long <- stats::reshape(
      tp_wide,
      direction = "long",
      varying   = names(tp_wide)[-1],
      v.names   = "sample_id",
      times     = names(tp_wide)[-1],
      timevar   = "timepoint",
      idvar     = "subject_id"
    )

    # remove the last variable, rename the first and reset the rownames
    tp_long   <- tp_long[, -ncol(tp_long)]
    names(tp_long)[names(tp_long) == "ind_id"] <- "subject_id"
    row.names(tp_long) <- NULL

    # filter the NAs out
    tp_long   <- tp_long[!is.na(tp_long$sample_id), ]

    # ---- reconcile comparisons --------------------------------------------
    if (!is.null(comparisons)) {
      # ---------- 1. reference sets -------------------------------------------
      valid_vals <- union(tp_long$timepoint,
                          unique(samples$group %||% character()))

      # ---------- 2. sanity-check comparisons ---------------------------------
      bad <- setdiff(
        unique(c(comparisons$group1, comparisons$group2)),
        valid_vals
        )

      .chk_cond(length(bad) > 0,
                sprintf("Comparisons refer to unknown group/timepoint: %s",
                        paste(bad, collapse = ", ")))

      # if group duplicates timepoint drop it
      if ("group" %in% names(samples) &&
          identical(samples$group,
                    tp_long$timepoint[match(samples$sample_id,
                                            tp_long$sample_id)])) {
        samples$group <- NULL
        comparisons$variable <- "timepoint"
      }
    }
    # ---------- 3. merge ------------------------------------------------------
    # add the time-point info if present
    samples <- merge(samples,
                     tp_long,
                     by = "sample_id"
    )

    timepoints <- tp_long
  }

  list(samples     = samples,
       comparisons = comparisons,
       timepoints  = timepoints)
}
