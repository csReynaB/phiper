# run heavy data prep only locally
# control with env var phiper_run_heavy_tests=1
# use phiper_prepare_dir from .renviron to locate the etl script

run_heavy <- Sys.getenv("PHIPER_RUN_HEAVY_TESTS", "0") == "1"

prep_dir <- Sys.getenv("PHIPER_PREPARE_DIR", unset = "")
if (!nzchar(prep_dir)) {
  .ph_warn(
    headline = "[tests] phiper_prepare_dir not set",
    step = "falling back to package-root/local-etl",
    bullets = c(
      "set PHIPER_PREPARE_DIR in .Renviron for a fixed absolute path",
      sprintf("cwd: %s", getwd())
    )
  )
  # fallback for ad hoc runs
  pkg_root <- tryCatch(
    rprojroot::find_package_root_file(),
    error = function(e) normalizePath(file.path(testthat::test_path("..")), mustWork = FALSE)
  )
  prep_dir <- file.path(pkg_root, "local-etl")
}

src <- normalizePath(file.path(prep_dir, "prepare_babies.R"), mustWork = FALSE)

if (run_heavy) {
  if (file.exists(src)) {
    .ph_with_timing(
      headline = "[tests] run local etl",
      step     = sprintf("sourcing: %s", src),
      expr     = source(src, local = new.env(parent = globalenv()))
    )
  } else {
    .ph_warn(
      headline = "[tests] etl script not found",
      step = "skipping heavy data preparation",
      bullets = c(
        sprintf("path tried: %s", src),
        sprintf("prep_dir: %s", prep_dir)
      )
    )
  }
} else {
  .ph_log_info(
    headline = "[tests] heavy tests off",
    bullets  = c("phiper_run_heavy_tests != 1", "or running on cran")
  )
}
