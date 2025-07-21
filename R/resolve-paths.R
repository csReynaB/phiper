#' @title Resolve legacy-import paths and perform fast-fail argument checks
#'
#' @description Combines explicit arguments with a YAML config (if given),
#'   expands every relative path to an absolute path (relative paths are
#'   evaluated against `dirname(config_yaml)` (!!!) when YAML is used, otherwise
#'   against the directory that contains the first supplied data matrix (!!!)),
#'   and returns a fully populated list of file locations and options ready for
#'   backend builders. Only cheap, load-blocking checks are done here:
#'
#' * `input_file` and `hit_file` must be supplied together or both omitted.
#'
#' * At least one matrix source (`exist_file`, `fold_change_file`, or the
#'   `input_file` + `hit_file` pair) must be present.
#'
#' * Deprecated `output_dir` triggers a soft warning.
#'
#'  All deeper table-content validation is deferred to `phip_data` class
#'  validation.
#'
#' @param exist_file,fold_change_file,input_file,hit_file,samples_file,timepoints_file,comparisons_file
#'  Character paths (relative or absolute) to the respective CSV/Parquet inputs.
#'  `NULL` means "not supplied".
#' @param extra_cols Character vector of extra metadata columns to keep; may be
#'   `NULL`.
#' @param output_dir Ignored (soft-deprecated).
#' @param backend Desired backend string (`"duckdb"`, `"arrow"`, `"memory"`);
#'   may come from YAML or explicit arg.
#' @param config_yaml Optional path to a YAML file whose keys mirror the
#'   function arguments; relative paths inside the YAML are resolved against the
#'   YAML’s own directory.
#'
#' @return A named list with absolute paths, `backend`, `extra_cols`, and
#'   `base_dir`; suitable for downstream helper functions.
#'
#' @keywords internal

.resolve_paths <- function(
    exist_file = NULL,
    fold_change_file = NULL,
    samples_file = NULL,
    input_file = NULL,
    hit_file = NULL,
    timepoints_file = NULL,
    extra_cols = NULL,
    comparisons_file = NULL,
    output_dir = NULL, # deprecated
    backend = NULL,
    config_yaml = NULL) {
  ## ------------------------------------------------------------------------ ##
  ## 1.  locate base directory & read yaml (if any provided)                  ##
  ## ------------------------------------------------------------------------ ##
  base_dir <- if (!is.null(config_yaml)) {
    fs::path_dir(fs::path_abs(config_yaml))
  } else {
    ## at least samples_file must be given when no YAML ----------------
    .chk_cond(
      is.null(samples_file),
      "When 'config_yaml' is absent you must supply 'samples_file'."
    )

    fs::path_dir(fs::path_abs(exist_file %||% samples_file))
  }

  yaml_cfg <- if (!is.null(config_yaml)) {
    ## validate the file extension
    .chk_path(config_yaml, "config_yaml", c("yml", "yaml"))

    ## read the yaml
    rlang::check_installed("yaml")
    yaml::read_yaml(config_yaml)
  } else {
    ## safe fallback
    NULL
  }

  ## ------------------------------------------------------------------------ ##
  ## 2.  helper to merge yaml + explicit args, validate & absolutise          ##
  ## ------------------------------------------------------------------------ ##
  fetch <- function(arg,
                    key,
                    validate = NULL,
                    optional = FALSE,
                    ...) {
    # safe fallback to NULL --> if both NULL, then %||% returns NULL
    val <- yaml_cfg[[key]] %||% arg # yaml first, then explicit

    # required argument have to be provided!
    .chk_cond(
      is.null(val) && !optional,
      sprintf("Missing required argument '%s' in YAML or call.", key)
    )

    # if the validator if .chk_path, expand the path to absolute
    if (!is.null(val) && identical(validate, .chk_path)) {
      if (!fs::is_absolute_path(val)) {
        val <- fs::path_abs(val, start = base_dir)
      }
    }

    # perform the custom validation if specified
    if (!is.null(val) && is.function(validate)) {
      validate(val, key, ...)
    }

    # return
    val
  }

  ## ------------------------------------------------------------------------ ##
  ## 3.  resolve every supported argument                                     ##
  ## ------------------------------------------------------------------------ ##
  cfg <- list(
    exist_file = fetch(exist_file,
      "exist_file",
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    fold_change_file = fetch(fold_change_file,
      "fold_change_file",
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    input_file = fetch(input_file,
      "input_file",
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    hit_file = fetch(hit_file,
      "hit_file",
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    samples_file = fetch(samples_file,
      "samples_file",
      .chk_path,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    timepoints_file = fetch(timepoints_file,
      "timepoints_file",
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    comparisons_file = fetch(comparisons_file,
      "comparisons_file",
      .chk_path,
      optional = TRUE,
      extension = c("csv", "parquet", "parq", "pq")
    ),
    extra_cols = fetch(extra_cols,
      "extra_cols",
      optional = TRUE
    ),
    output_dir = fetch(output_dir,
      "output_dir",
      optional = TRUE
    ),
    backend = fetch(
      backend,
      "backend",
      chk::chk_string
    ),
    base_dir = base_dir # for downstream helpers
  )

  ## ------------------------------------------------------------------------ ##
  ## 4.  fast-fail rules that really must hold before heavy work              ##
  ## ------------------------------------------------------------------------ ##
  #  rule 1: input_file and hit_file must be provided together ----------
  .chk_cond(
    xor(is.null(cfg$input_file), is.null(cfg$hit_file)),
    "Arguments 'input_file' and 'hit_file' must be provided together."
  )

  #  rule 2: at least one matrix source present -------------------------
  all_null <- with(cfg, is.null(exist_file)  &&
                     is.null(fold_change_file) &&
                     is.null(input_file) &&
                     is.null(hit_file))
  .chk_cond(all_null,
            paste0("Supply at least one of:\n",
                   "  * 'exist_file'\n",
                   "  * 'fold_change_file'\n",
                   "  * both 'input_file' and 'hit_file'"))

  #  deprecation notice -------------------------------------------------
  .chk_cond(!is.null(cfg$output_dir),
    error = FALSE,
    "'output_dir' is deprecated and will be ignored."
  )

  # validate the tables itself --> the logic has been moved entirely to the
  # phip_data class validator

  cfg
}
