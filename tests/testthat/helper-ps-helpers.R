# tests/testthat/helper-ps-helpers.R
# helpers to access heavy objects exposed via options by the etl

get_heavy_ps <- function() getOption("phiper.test.ps", NULL)
get_heavy_meta_kids <- function() getOption("phiper.test.meta_kids", NULL)

skip_if_no_heavy <- function() {
  if (is.null(get_heavy_ps())) {
    testthat::skip("heavy data not available enable PHIPER_RUN_HEAVY_TESTS=1
                   and ensure prepare_babies.R sets options")
  }
  invisible(TRUE)
}
